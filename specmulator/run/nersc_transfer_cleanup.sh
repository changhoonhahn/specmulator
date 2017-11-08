#!bin/bash 
source_dir=/Users/chang/projects/specmulator/dat
dest_dir=/global/homes/c/chahah/projects/specmulator/dat

echo "password"
read -s pwd
echo "overwrite? y/n"
read owrite

for mneut in 0.15eV  
do 
    # check that the neutrino directory exists
    mneutdir=$dest_dir/$mneut/
    if sshpass -p $pwd ssh edison '[ ! -d '$mneutdir' ]'; then
        sshpass -p $pwd ssh edison "mkdir "$mneutdir
    fi 

    for real in 2 89 90; do 
        groupfile_source=$source_dir/$mneut/$real/groups_004/group_ids_004.0
        if [ -f $groupfile_source ]; then 
            # check that the realization directory exists
            realdir=$dest_dir/$mneut/$real/
            if sshpass -p $pwd ssh edison '[ ! -d '$realdir' ]'; then 
                sshpass -p $pwd ssh edison "mkdir "$realdir
            fi 
            # check that the snapshot directory exists
            snapdir=$dest_dir/$mneut/$real/snapdir_004/
            if sshpass -p $pwd ssh edison '[ ! -d '$snapdir' ]'; then 
                sshpass -p $pwd ssh edison "mkdir "$snapdir
            fi 
            # check whether snapshot file exists
            snapfile=$snapdir'snap_004.0'
            if sshpass -p $pwd ssh edison '[ ! -f '$snapfile' ]'; then
                sshpass -p $pwd scp $source_dir/$mneut/$real/snapdir_004/snap_004.0 edison:$snapdir
            else 
                if [ $owrite = "y" ]; then 
                    echo "overwriting "$snapfile
                    sshpass -p $pwd ssh edison "rm "$snapfile
                    sshpass -p $pwd scp $source_dir/$mneut/$real/snapdir_004/snap_004.0 edison:$snapdir
                else
                    echo $snapfile" exists" 
                fi 
            fi 
            # check whether group file exists
            groupfile=$dest_dir/$mneut/$real/groups_004/group_ids_004.0
            groupdir=$dest_dir/$mneut/$real/groups_004/
            if sshpass -p $pwd ssh edison '[ ! -f '$groupfile' ]'; then
                sshpass -p $pwd scp -r $source_dir/$mneut/$real/groups_004/ edison:$groupdir
            else 
                if [ $owrite = "y" ]; then 
                    echo "overwriting "$groupfile 
                    echo sshpass -p $pwd scp -r $source_dir/$mneut/$real/groups_004/ edison:$groupdir
                    sshpass -p $pwd scp -r $source_dir/$mneut/$real/groups_004/ edison:$groupdir
                else
                    echo $groupfile" exists" 
                fi 
            fi 
        else
            echo $groupfile_source" does not exist" 
        fi 
    done
done 
