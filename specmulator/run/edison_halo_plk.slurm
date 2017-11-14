#!/bin/bash
#SBATCH -p debug 
#SBATCH -n 1
#SBATCH -t 00:10:00 
#SBATCH -J plk_halo_0.6eV_44
#SBATCH -o plk_halo_0.6eV_44.o

source /usr/common/contrib/bccp/conda-activate.sh 3.6 

codedir=/global/homes/c/chahah/projects/specmulator/specmulator

for mneut in 0.6; do
    for real in 44; do #{1..100}; do 
        srun -n 1 -c 1 python $codedir/run/run_halo_catalog.py plk $mneut $real 4 real
        srun -n 1 -c 1 python $codedir/run/run_halo_catalog.py plk $mneut $real 4 z 
    done
done 
