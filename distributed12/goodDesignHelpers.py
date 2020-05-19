import numpy as np
import autograd.numpy as np
from autograd import grad, jacobian
from numpy import heaviside
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from collections import Counter
from functools import partial
import json
from math import pow, exp
from scipy.optimize import fsolve, root
from functools import partial


# heavidside function
def Heaviside(x,theta):
    return 1 if x-theta >=0 else 0

# hill function
def Hill(x,theta,v):
    return 1/((theta/x)**v+ 1)

# function that returns dz/dt
def model(z,t,H,theta1,theta2):
    dm1dt = 1 - mu*z[0] - H(z[3],theta1)*z[0] + H(z[2],theta1)*z[1]
    dm2dt = 1 - mu*z[1] + H(z[3],theta1)*z[0] - H(z[2],theta1)*z[1]
    dp1dt = (beta+H(z[2],theta2))*z[0] - pi*z[2] + D*(z[3]-z[2])
    dp2dt = (beta+H(z[3],theta2))*z[1] - pi*z[3] + D*(z[2]-z[3])
    dzdt = [dm1dt,dm2dt,dp1dt,dp2dt]
    return dzdt


def isInRegion(p,params):
    p1 ,p2 = p
    theta1,theta2, mu,beta,pi,d = params
    
    first = (2*beta/(mu*pi) <= (p1+p2) <= 2*(beta+1)/(mu*pi))
    second = (abs(p1-p2) <=(2+2*beta)/(mu*(pi+d)))
    
    if first and second:
        return True
    else:
        False
        
        
def isAttractorOverall(p, parameters,H1,H2):
    theta1, theta2, mu, beta, pi, d = parameters
    
    def computeFixedPoint(p):
        p1,p2 = p
        m1 = (1+2/mu*H1(p1,theta1))/(mu+H1(p2,theta1)+H1(p1,theta1))
        m2 = 2/mu - m1
        return [m1,m2,p1,p2]
    
    fp = computeFixedPoint(p)
    
    def M1(fp):
        m1,m2,p1,p2 = fp
        return 1-mu*m1 - H1(p2,theta1)*m1 + H1(p1,theta1)*m2

    def M2(fp):
        m1,m2,p1,p2 = fp
        return 1-mu*m2 + H1(p2,theta1)*m1 - H1(p1,theta1)*m2

    def P1(fp):
        m1,m2,p1,p2 = fp
        return (beta+ H2(p1,theta2))*m1 - pi*p1 + d*(p2-p1)

    def P2(fp):
        m1,m2,p1,p2 = fp
        return (beta+ H2(p2,theta2))*m2 - pi*p2 + d*(p1-p2)
    
    matrix = [grad(M1)(np.array(fp)),grad(M2)(np.array(fp)),grad(P1)(np.array(fp)),grad(P2)(np.array(fp))]
    vals = np.linalg.eigvals(matrix)
    
    for v in vals:
        if v>=0:
            return False
    return True

def regionSample(thetaHigh,thetaLow,bound):
    lower = list(np.linspace(0,thetaLow,11)[1:])
    middle = list(np.linspace(thetaLow,thetaHigh,11)[1:])
    bigger = list(np.linspace(thetaHigh,max(bound,2*thetaHigh),10))
    
    return lower+middle+bigger


# algebraic expression of fixed point for p1 and p2
# are we going to check the eigenvalues is negative or not
def findFP(p,params,H1,H2):
    
    p1 ,p2 = p
    theta1,theta2, mu,beta,pi,d = params
    
    v1 = (beta+H2(p1,theta2))*(1+2/mu*H1(p1,theta1))/(mu+H1(p1,theta1)+H1(p2,theta1)) - (pi+d)*p1 + d*p2
    
    v2 = (beta+H2(p2,theta2))*2/mu - (beta+H2(p2,theta2))*(1+2/mu*H1(p1,theta1))/(mu+H1(p2,theta1)+H1(p1,theta1))\
    +d*p1 - (pi+d)*p2
    
    return v1, v2


# portraits is a dictionary of number fixed point in the corresponding 9 regions
def generatePortraits(p,parameters,portrait):
    theta1, theta2, mu, beta, pi, d = parameters
    p1,p2 = p
    
    y, x = 0, 0

    if p1 > theta1:
        if p1 < theta2:
            x = 1
        else:
            x = 2

    if p2 > theta1:
        if p2 < theta2:
            y = 1
        else:
            y = 2

    val = 3*y + x
    
    if val not in portrait:
        portrait[val] = set()
        
    roundP = tuple(map(lambda x: round(x,6),p))
    
    portrait[val].add(tuple(roundP))
    
def encode(portrait):
    ret = []
    for i in [3,6,7,0,4,8]:
        if i in portrait:
            ret.append(len(portrait[i]))
        else:
            ret.append(0)
    
    return ret