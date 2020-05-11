%% 2 cell orb mRNA/protein distribution
Delta = 1; %symmetric mRNA transport

%parameter set is given (theta_1, theta_2, mu, beta, pi, D)
param_set =  [1.8, 0.67, 0.275, 0.45, 2.95, 0.8];
theta1 = param_set(1);
theta2 = param_set(2);
mu = param_set(3);
beta = param_set(4);
pi = param_set(5);
D = param_set(6);

%fxns for theta_2 < theta_1 case to check which region the param set is in
r0 = beta/(mu*pi);
r1 = (beta*(pi+2*D)+D)/(mu*pi*(pi+2*D));
r2 = (beta*(pi+2*D)+(pi+D))/(mu*pi*(pi+2*D));
r3 = (beta*mu*D+(beta+1)*(mu+2)*(pi+D))/(mu*(1+mu)*pi*(pi+2*D));
r4 = (beta*mu*(pi+D)+(beta+1)*(mu+2)*D)/(mu*(1+mu)*pi*(pi+2*D));
r5 = (beta+1)/(mu*pi);
r6 = ((beta+1)*(mu*(pi+2*D)+2*(pi+D)))/(mu*(1+mu)*pi*(pi+2*D));
r7 = ((beta+1)*(mu*(pi+2*D)+2*D))/(mu*(1+mu)*pi*(pi+2*D));

%confirm given parameter set lies in the region 18:1 and look at dynamics
check_18_1 = r5 < theta1 & theta1 < r3;
check_18_2 = max(r0,r4) < theta2 & theta2 < min(r1,r7);

%% test random ICs to get example trajectories
m1_init = 5*rand(1);
m2_init = 5*rand(1);
p1_init = [0 0.5 0.75 1 1.25 1.5 2 3.5 5];
p2_init = [0 0.5 0.75 1 1.25 1.5 2 3.5 5];

%can plot these lines to see param regionsin phase plane
y = linspace(0,10,100);
theta1y = theta1*ones(length(y),1);
theta2y = theta2*ones(length(y),1);

%pick finite Hill coefficients
n1 = 20;
n2 = 20;

figure;
titlestring = strcat(['$(\theta_1,\theta_2,\mu,\beta,\pi,D, n, \nu) = (',num2str(theta1),',',num2str(theta2),',',num2str(mu),',',num2str(beta),',',num2str(pi),',',num2str(D),',',num2str(n1),',',num2str(n2),')$']);
title1 = strcat(['$n_1,n_2 = ',num2str(n1),',',num2str(n2),'$']);
hold on
for i = 1:length(p2_init)
    for j = 1:length(p1_init)
        x0 = [m1_init p1_init(j) m2_init p2_init(i)];
        
        [t,x] = ode23s(@(t,x) [(1-mu*x(1)-x(4).^n1/(theta1.^n1+x(4).^n1)*x(1)+Delta*x(2).^n1/(theta1.^n1+x(2).^n1)*x(3));
            ((beta+x(2).^n2/(x(2).^n2+theta2.^n2))*x(1)-pi*x(2)+D*(x(4)-x(2)));
            (1-mu*x(3)+x(4).^n1/(theta1.^n1+x(4).^n1)*x(1)-Delta*x(2).^n1/(theta1.^n1+x(2).^n1)*x(3));
            ((beta+x(4).^n2/(x(4).^n2+theta2.^n2))*x(3)-pi*x(4)+D*(x(2)-x(4)))], [0 10000], x0);
        
        plot(x(:,2),x(:,4),'k-','LineWidth',1)
        plot(x(end,2),x(end,4),'rx','MarkerSize',24,'LineWidth',4)
        
    end
end
plot(linspace(0,5,1000),linspace(0,5,1000),'k--','LineWidth',4)
xlabel('$P_1$','interpreter','latex','FontSize',24)
ylabel('$P_2$','interpreter','latex','FontSize',24)
% title(titlestring,'interpreter','latex')
% grid on
box on
axis square
axis([0 5 0 5])

%% Check stability case by case
%define steady state given param set through gradient descent
steady_state  = [];

theta1 = param_set(1);
theta2 = param_set(2);
mu = param_set(3);
beta = param_set(4);
pi = param_set(5);
D = param_set(6);

m1 = steady_state(1);
m2 = steady_state(2);
p1 = steady_state(3);
p2 = steady_state(4);

n1 = 20;
n2 = 20;

f1 = p1.^n1./(p1.^n1+theta1.^n1);
f2 = p2.^n1./(p2.^n1+theta1.^n1);
g1 = p1.^n2./(p1.^n2+theta2.^n2);
g2 = p2.^n2./(p2.^n2+theta2.^n2);

fp1 = n1.*f1.*(1-f1)./p1;
fp2 = n1.*f2.*(1-f2)./p2;
gp1 = n2.*g1.*(1-g1)./p1;
gp2 = n2.*g2.*(1-g2)./p2;

J = [-mu-f2 f1 fp1*m2 -fp2*m1;
    f2 -mu-f1 -fp1*m2 fp2*m1;
    beta+g1 0 gp1*m1-pi-D D;
    0 beta+g2 D gp2*m2-pi-D];

A = [(1-mu*m1-p2.^n1/(theta1.^n1+p2.^n1)*m1+p1.^n1/(theta1.^n1+p1.^n1)*m2);
     (1-mu*m2+p2.^n1/(theta1.^n1+p2.^n1)*m1-p1.^n1/(theta1.^n1+p1.^n1)*m2);
     ((beta+p1.^n2/(p1.^n2+theta2.^n2))*m1-pi*p1+D*(p2-p1));
     ((beta+p2.^n2/(p2.^n2+theta2.^n2))*m2-pi*p2+D*(p1-p2))];

eig(J) %look at eigenvalues to check stability
disp(sum(abs(A))); %confirm this is a steady state