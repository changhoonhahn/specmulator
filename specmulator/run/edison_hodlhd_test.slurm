#!/bin/bash
#SBATCH -p debug 
#SBATCH -n 1
#SBATCH -t 00:30:00 
#SBATCH -J hodlhd_test
#SBATCH -o hodlhd_test.o

source /usr/common/contrib/bccp/conda-activate.sh 3.6

mneut=0.0
nreal=1
nzbin=4

for rseed in 1; do 
    for i_hod in 0; do 
        srun -n 1 -c 1 python /global/homes/c/chahah/projects/specmulator/specmulator/run/run_hodlhd_catalog.py catalog 40 $mneut $nreal $nzbin $rseed $i_hod
        srun -n 1 -c 1 python /global/homes/c/chahah/projects/specmulator/specmulator/run/run_hodlhd_catalog.py observable 40 $mneut $nreal $nzbin $rseed $i_hod plk real 
        srun -n 1 -c 1 python /global/homes/c/chahah/projects/specmulator/specmulator/run/run_hodlhd_catalog.py observable 40 $mneut $nreal $nzbin $rseed $i_hod plk z 
    done 
done 
