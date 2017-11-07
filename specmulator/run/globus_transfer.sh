#!bin/bash/
source_ep='7aec0828-e9fb-11e5-97d6-22000b9da45e'
dest_ep='1c93ed5e-b8da-11e7-b143-22000a92523b'

sourdir=/mnt/ceph/users/fvillaescusa/Neutrino_simulations/Sims_Dec16_2
destdir=/Users/chang/projects/specmulator/dat

for mneut in 0.9eV # 0.10eV 0.15eV 0.6eV
do 
    for real in {1..2}
    do 
        # check whether the realization directory exists 
        realdir=$destdir/$mneut/$real
        [ -d $realdir ] || mkdir $realdir
        # check whether the realization directory exists 
        groupdir=$destdir/$mneut/$real/groups_004
        [ -d $groupdir ] || mkdir $groupdir 
        # check whether the snapshot directory exists 
        snapdir=$destdir/$mneut/$real/snapdir_004
        [ -d $snapdir ] || mkdir $snapdir 

        snapfile=$destdir/$mneut/$real/snapdir_004/snap_004.0
        if [ ! -f $snapfile ];
        then 
            globus transfer $source_ep:$sourdir/$mneut/$real/snapdir_004/snap_004.0 \
                $dest_ep:$destdir/$mneut/$real/snapdir_004/snap_004.0
        else
            echo $snapfile' already exists'
        fi 
        groupfile=$destdir/$mneut/$real/groups_004/group_ids_004.0
        if [ ! -f $groupfile ]; 
        then 
            globus transfer $source_ep:$sourdir/$mneut/$real/groups_004/ \
                $dest_ep:$destdir/$mneut/$real/groups_004 --recursive
            sleep 120s 
        else
            echo $groupfile' already exists'
        fi 
    done 
done 
