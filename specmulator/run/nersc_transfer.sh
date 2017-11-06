#!bin/bash 
source_dir=/Users/chang/projects/specmulator/dat
dest_dir=/global/homes/c/chahah/projects/specmulator/dat

echo "password"
read -s pwd

for mneut in 0.0eV 0.06eV 0.10eV 0.15eV 0.6eV
do 
    for real in {3..10} 
    do 
        # check that the neutrino directory exists
        mneutdir=$dest_dir/$mneut/
        if sshpass -p $pwd ssh edison '[ ! -d '$snapdir' ]';
        then
            sshpass -p $pwd ssh edison "mkdir "$mneutdir
        fi 
        realdir_source=$source_dir/$mneut/$real/
        if 
        then 
            # check that the realization directory exists
            realdir=$dest_dir/$mneut/$real/
            if sshpass -p $pwd ssh edison '[ ! -d '$realdir' ]';
            then 
                sshpass -p $pwd ssh edison "mkdir "$realdir
            fi 
            # check that the snapshot directory exists
            snapdir=$dest_dir/$mneut/$real/snapdir_004/
            if sshpass -p $pwd ssh edison '[ ! -d '$snapdir' ]';
            then 
                sshpass -p $pwd ssh edison "mkdir "$snapdir
            fi 
            # check whether snapshot file exists
            snapfile=$snapdir'snap_004.0'
            if sshpass -p $pwd ssh edison '[ ! -f '$snapfile' ]';
            then
                sshpass -p $pwd scp $source_dir/$mneut/$real/snapdir_004/snap_004.0 edison:$snapdir
            fi 
            # check whether group file exists
            groupfile=$dest_dir/$mneut/$real/groups_004/group_ids_004.0
            if sshpass -p $pwd ssh edison '[ ! -f '$groupfile' ]';
            then
                groupdir=$dest_dir/$mneut/$real/groups_004/
                sshpass -p $pwd scp -r $source_dir/$mneut/$real/groups_004/ edison:$groupdir
            fi 
        fi 
    done
done 
