#!/bin/bash

# Active comments for SLURM
#SBATCH -n 1                  # One task
#SBATCH -c 1                 # One cpu per task
#SBATCH -N 1                   # Minimum one node
#SBATCH -t 0-12:00:00            # Runtime in D-HH:MM
#SBATCH -p mischaik_1 # Partition to submit to
# #SBATCH --mem-per-cpu=4000   # Memory pool for all cores (see also --mem-per-cpu)

# Optional arguments (uncomment them to use)
#SBATCH --output=testslurm%N.%j.out    # Output file
#SBATCH --error=slurm_script.%N.%j.err     # Error output file
#SBATCH --mail-user=chamberlian1990@gmail.com  # User e-mail
#SBATCH --mail-type=FAIL         # When to send e-mail
#SBATCH --mem=8G # 16GB of memory
# #SBATCH --array=1-3

srun python gdAmarel.py  1
srun python gdAmarel.py  2
srun python gdAmarel.py  3
srun python gdAmarel.py  4
srun python gdAmarel.py  5
srun python gdAmarel.py  6
srun python gdAmarel.py  7
srun python gdAmarel.py  8
srun python gdAmarel.py  9
srun python gdAmarel.py  10
srun python gdAmarel.py  11
srun python gdAmarel.py  12
srun python gdAmarel.py  13
srun python gdAmarel.py  14
srun python gdAmarel.py  15
srun python gdAmarel.py  16
srun python gdAmarel.py  17
srun python gdAmarel.py  18
srun python gdAmarel.py  19
srun python gdAmarel.py  20
srun python gdAmarel.py  21
srun python gdAmarel.py  22
srun python gdAmarel.py  23
srun python gdAmarel.py  24
srun python gdAmarel.py  25
srun python gdAmarel.py  26
srun python gdAmarel.py  27
srun python gdAmarel.py  28
srun python gdAmarel.py  29
srun python gdAmarel.py  30
srun python gdAmarel.py  31
srun python gdAmarel.py  32
srun python gdAmarel.py  33
srun python gdAmarel.py  34
srun python gdAmarel.py  35
srun python gdAmarel.py  36
srun python gdAmarel.py  37
srun python gdAmarel.py  38
srun python gdAmarel.py  39
