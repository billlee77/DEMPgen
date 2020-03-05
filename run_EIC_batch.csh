#!/bin/csh

set FileNum=$1
set NumEvents=$2
echo "Running target polarisation up, FF setting for file $FileNum with $NumEvents events per file"
sed -i 's/"file_name" \:.*/"file_name" \: "DEMPGen_'$NumEvents'_'$FileNum'",/' Config_EIC.json
sed -i 's/"n_events" \:.*/"n_events" \: '$NumEvents',/' Config_EIC.json 
cd data/
./../build/DEMPgen ../Config_EIC.json
sleep 5
mv "LundFiles/eic_input_DEMPGen_"$NumEvents"_"$FileNum".dat" "LundFiles/eic_DEMPGen_"$NumEvents"_"$FileNum".lund"
