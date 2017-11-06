#!/bin/bash
#SBATCH -p debug 
#SBATCH -n 1
#SBATCH -t 00:10:00 
#SBATCH -J pre3PCF
#SBATCH -o pre3PCF.o

source /usr/common/contrib/bccp/conda-activate.sh 3.6 

srun -n 1 -c 1 python /global/homes/c/chahah/projects/specmulator/specmulator/run/run_halo_catalog.py 3pcf 0.0 1 4 real
srun -n 1 -c 1 python /global/homes/c/chahah/projects/specmulator/specmulator/run/run_halo_catalog.py 3pcf 0.0 1 4 z 
