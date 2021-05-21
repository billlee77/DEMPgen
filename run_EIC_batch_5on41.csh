#!/bin/csh

set FileNum=$1
set NumEvents=$2
set RandomSeed=$3
echo "Running target polarisation up, FF setting for file $FileNum with $NumEvents events per file using random seed $RandomSeed"
cp Config_EIC.json Config_EIC_5on41_$FileNum.json
sed -i 's/"file_name" \:.*/"file_name" \: "DEMPGen_5on41_'$NumEvents'_'$FileNum'",/' Config_EICC_5on41_$FileNum.json
sed -i 's/"n_events" \:.*/"n_events" \: '$NumEvents',/' Config_EICC_5on41_$FileNum.json
sed -i 's/"generator_seed"\:.*/"generator_seed" \: '$RandomSeed',/' Config_EICC_5on41_$FileNum.json
cd data/
./../build/DEMPgen ../Config_EICC_5on41_$FileNum.json
sleep 5
mv "LundFiles/eic_input_DEMPGen_5on41_"$NumEvents"_"$FileNum".dat" "LundFiles/eic_DEMPGen_5on41_"$NumEvents"_"$FileNum".lund"
rm -rf ../Config_EICC_5on41_$FileNum.json
