#!/bin/bash
#SBATCH -p debug 
#SBATCH -n 1
#SBATCH -t 00:05:00 
#SBATCH -J pre3PCF
#SBATCH -o pre3PCF.o

source /usr/common/contrib/bccp/conda-activate.sh 3.6 
code_dir=/global/homes/c/chahah/projects/specmulator/specmulator

for mneut in 0.6; do
    srun -n 1 -c 1 python $code_dir/run/run_halo_catalog.py 3pcf $mneut 44 44 4 real 1
    srun -n 1 -c 1 python $code_dir/run/run_halo_catalog.py 3pcf $mneut 44 44 4 z 1 
done
