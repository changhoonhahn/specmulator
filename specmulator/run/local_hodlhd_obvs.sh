#!/bin/bash

mneut=0.0
nreal=1
nzbin=4
rseed=1

for i_hod in {0..16}; do 
    python run/run_hodlhd_catalog.py observable $mneut $nreal $nzbin $rseed $i_hod plk real
    python run/run_hodlhd_catalog.py observable $mneut $nreal $nzbin $rseed $i_hod plk z 
done 
