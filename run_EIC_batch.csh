#!/bin/csh

set FileNum=$1
set NumEvents=$2
set RandomSeed=$3
echo "Running target polarisation up, FF setting for file $FileNum with $NumEvents events per file using random seed $RandomSeed"
cp Config_EIC.json Config_EIC_$FileNum.json
sed -i 's/"file_name" \:.*/"file_name" \: "DEMPGen_5on100_'$NumEvents'_'$FileNum'",/' Config_EIC_$FileNum.json
sed -i 's/"n_events" \:.*/"n_events" \: '$NumEvents',/' Config_EIC_$FileNum.json
sed -i 's/"generator_seed"\:.*/"generator_seed" \: '$RandomSeed',/' Config_EIC_$FileNum.json
cd data/
./../build/DEMPgen ../Config_EIC_$FileNum.json
sleep 5
mv "LundFiles/eic_input_DEMPGen_5on100_"$NumEvents"_"$FileNum".dat" "LundFiles/eic_DEMPGen_5on100_"$NumEvents"_"$FileNum".lund"
rm -rf ../Config_EIC_$FileNum.json
