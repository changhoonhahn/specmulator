#!/bin/bash
#SBATCH -p regular 
#SBATCH -n 1
#SBATCH -t 05:00:00 
#SBATCH -J hodlhd
#SBATCH -o hodlhd.o

source /usr/common/contrib/bccp/conda-activate.sh 3.6 

mneut=0.0
nreal=1
nzbin=4
rseed=1

for i_hod in {0..99}; do 
    srun -n 1 -c 1 python /global/homes/c/chahah/projects/specmulator/specmulator/run/run_hodlhd_catalog.py catalog 100 $mneut $nreal $nzbin $rseed $i_hod
    srun -n 1 -c 1 python /global/homes/c/chahah/projects/specmulator/specmulator/run/run_hodlhd_catalog.py observable 100 $mneut $nreal $nzbin $rseed $i_hod plk real 
    srun -n 1 -c 1 python /global/homes/c/chahah/projects/specmulator/specmulator/run/run_hodlhd_catalog.py observable 100 $mneut $nreal $nzbin $rseed $i_hod plk z 
done 
