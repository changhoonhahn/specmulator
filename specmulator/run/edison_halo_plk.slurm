#!/bin/bash
#SBATCH -p debug 
#SBATCH -n 1
#SBATCH -t 00:30:00 
#SBATCH -J plk_halo 
#SBATCH -o plk_halo.o

source /usr/common/contrib/bccp/conda-activate.sh 3.6 

codedir=/global/homes/c/chahah/projects/specmulator/specmulator

for mneut in 0.0; do
    for real in {2..10}; do 
        srun -n 1 -c 1 python $codedir/run/run_halo_catalog.py plk $mneut $real 4 real
        srun -n 1 -c 1 python $codedir/run/run_halo_catalog.py plk $mneut $real 4 z 
    done
done 
