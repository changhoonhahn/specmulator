#!bin/bash 
edison_dir=/scratch2/scratchdirs/chahah/specmulator/lhd/0.0eV_1_z4_40samples/
local_dir=/Volumes/chang_eHDD/projects/specmulator/lhd/0.0eV_1_z4_40samples/

echo "password"
read -s pwd

for rseed in 1; do
    for i_hod in {0..39}; do 
        HOD_dir="HODmdu_seed"$rseed"_"$i_hod

        # check whether the directory exists
        if [ ! -d '$local_dir$HOD_dir' ]; then
            mkdir $local_dir$HOD_dir
        fi 
        sshpass -p $pwd scp edison:$edison_dir$HOD_dir/pk* $local_dir$HOD_dir/
    done 
done 
