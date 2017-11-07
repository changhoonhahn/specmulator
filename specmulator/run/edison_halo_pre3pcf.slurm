#!/bin/bash
#SBATCH -p regular 
#SBATCH -n 1
#SBATCH -t 00:30:00 
#SBATCH -J pre3PCF
#SBATCH -o pre3PCF.o

source /usr/common/contrib/bccp/conda-activate.sh 3.6 
code_dir=/global/homes/c/chahah/projects/specmulator/specmulator

for mneut in 0.0; do
    srun -n 1 -c 1 python $code_dir/run/run_halo_catalog.py 3pcf $mneut 12 100 4 real 10 
    srun -n 1 -c 1 python $code_dir/run/run_halo_catalog.py 3pcf $mneut 12 100 4 z 10 
done
