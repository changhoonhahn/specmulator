#!/bin/bash
#SBATCH -p regular 
#SBATCH -n 1
#SBATCH -t 00:30:00 
#SBATCH -J hodlhd
#SBATCH -o hodlhd.o

source /usr/common/contrib/bccp/conda-activate.sh 3.6 
#module load python/3.6-anaconda-4.4
#module load python/2.7-anaconda

#srun -n 1 -c 1 python /global/homes/c/chahah/projects/specmulator/specmulator/run/run_hodlhd_catalog.py lhd 
srun -n 1 -c 1 python /global/homes/c/chahah/projects/specmulator/specmulator/run/run_hodlhd_catalog.py catalog 0.0 1 4 1
