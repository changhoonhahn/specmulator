#!/bin/bash
#SBATCH -p debug 
#SBATCH -n 1
#SBATCH -t 00:15:00 
#SBATCH -J plk_halo 
#SBATCH -o plk_halo.o

source /usr/common/contrib/bccp/conda-activate.sh 3.6 

codedir=/global/homes/c/chahah/projects/specmulator/specmulator

for mneut in 0.0; do
    for real in {4..10}; do 
        now=$(date +"%T")
        echo $now
        srun -n 1 -c 1 python $codedir/run/run_halo_catalog.py plk $mneut $real 4 real
        now=$(date +"%T")
        echo $now
        srun -n 1 -c 1 python $codedir/run/run_halo_catalog.py plk $mneut $real 4 z 
        now=$(date +"%T")
        echo $now
    done
done 
