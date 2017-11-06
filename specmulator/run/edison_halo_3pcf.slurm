#!/bin/bash
#SBATCH -p debug 
#SBATCH -n 1
#SBATCH -t 00:20:00 
#SBATCH -J 3PCF_halo
#SBATCH -o 3PCF_halo.o

nside=80
dest_dir=/global/homes/c/chahah/projects/specmulator/dat

for mneut in 0.0eV; do
    for real in {2..10}; do 
        srun -n 1 -c 16 /global/homes/c/chahah/projects/specmulator/specmulator/zslepian/grid_multipoles \
            -box 1000. -scale 1 -nside $nside -in $dest_dir/$mneut/$real/groups.nzbin4.rspace.dat \
                > $dest_dir/$mneut/"3pcf_groups."$mneut"."$real".nzbin4.rspace.dat"
        srun -n 1 -c 16 /global/homes/c/chahah/projects/specmulator/specmulator/zslepian/grid_multipoles \
            -box 1000. -scale 1 -nside $nside -in $dest_dir/$mneut/$real/groups.nzbin4.zspace.dat \
                > $dest_dir/$mneut/"3pcf_groups."$mneut"."$real".nzbin4.zspace.dat"
    done 
done
