#!/bin/bash -l 
#SBATCH -p regular 
#SBATCH -N 1 
#SBATCH -t 00:30:00 
#SBATCH -J neuthalos 
#SBATCH -o neuthalos.o%j 
#SBATCH -L SCRATCH,project 

source /usr/common/contrib/bccp/conda-activate.sh 2.7 

srun -n 1 python /global/homes/c/chahah/projects/specmulator/specmulator/run/neutsims.py 0.0 2 4
