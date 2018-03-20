#!/bin/bash

mneut=0.0
nreal=1
nzbin=4

python run/run_hodlhd_catalog.py lhd 40

for rseed in 1; do #{1..10}; do 
    # overwrite LHD file 
    python run/run_hodlhd_catalog.py lhd 40
    for i_hod in {0..39}; do 
        python run/run_hodlhd_catalog.py catalog 40 $mneut $nreal $nzbin $rseed $i_hod
        python run/run_hodlhd_catalog.py observable 40 $mneut $nreal $nzbin $rseed $i_hod plk real
        python run/run_hodlhd_catalog.py observable 40 $mneut $nreal $nzbin $rseed $i_hod plk z 
    done 
done 
