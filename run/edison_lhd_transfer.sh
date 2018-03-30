#!bin/bash 
#edison_dir=/scratch2/scratchdirs/chahah/specmulator/lhd/0.0eV_1_z4_40samples/
#local_dir=/Volumes/chang_eHDD/projects/specmulator/lhd/0.0eV_1_z4_40samples/
edison_dir=/scratch2/scratchdirs/chahah/specmulator/lhd/onlyHOD/test_20_sinha2017prior_narrow/
local_dir=/Volumes/chang_eHDD/projects/specmulator/lhd/onlyHOD/test_20_sinha2017prior_narrow/

echo "password"
read -s pwd

for rseed in {1..10}; do
    for i_hod in {0..19}; do 
        #HOD_dir="HODmdu_seed"$rseed"_"$i_hod
        HOD_dir="testHOD_seed"$rseed"_"$i_hod
        echo $local_dir$HOD_dir
        # check whether the directory exists
        if [ ! -d '$local_dir$HOD_dir' ]; then
            mkdir -p $local_dir$HOD_dir
        fi 
        sshpass -p $pwd scp edison:$edison_dir$HOD_dir/pk* $local_dir$HOD_dir/
    done 
done 
