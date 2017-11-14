#!/bin/bash/

mneut=0.0
nreal=1
nzbin=4

for rseed in {1..100}; do 
    python run/run_hodlhd_catalog.py catalog 0 $mneut $nreal $nzbin $rseed fid 
    python run/run_hodlhd_catalog.py observable 0 $mneut $nreal $nzbin $rseed fid plk real
    python run/run_hodlhd_catalog.py observable 0 $mneut $nreal $nzbin $rseed fid plk z 
done 
