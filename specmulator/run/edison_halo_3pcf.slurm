#!/bin/bash
#SBATCH -p regular 
#SBATCH -n 1
#SBATCH -t 00:10:00 
#SBATCH -J 3PCF_halo_0.10eV
#SBATCH -o 3PCF_halo_0.10eV.o

nbin=20
nside=20
dest_dir=/global/homes/c/chahah/projects/specmulator/dat

for mneut in 0.10eV; do
    for real in 100; do #{1..100}; do 
        out3pcf=$dest_dir/$mneut/"3pcf.groups."$mneut"."$real".nzbin4.nside"$nside".nbin"$nbin".rspace.dat"
        srun -n 1 -c 16 /global/homes/c/chahah/projects/specmulator/specmulator/zslepian/grid_multipoles_nbin20 \
            -box 1000. -scale 1 -nside $nside -in $dest_dir/$mneut/$real/groups.nzbin4.rspace.dat > $out3pcf
        out3pcf_z=$dest_dir/$mneut/"3pcf.groups."$mneut"."$real".nzbin4.nside"$nside".nbin"$nbin".zspace.dat"
        srun -n 1 -c 16 /global/homes/c/chahah/projects/specmulator/specmulator/zslepian/grid_multipoles_nbin20 \
            -box 1000. -scale 1 -nside $nside -in $dest_dir/$mneut/$real/groups.nzbin4.zspace.dat > $out3pcf_z
    done 
done
