#!/bin/bash 

folder=RESLT_example 
logfile=$folder"/oomph-out.txt"
parameters="--reynolds 44.795569 --gravity 0 --capilliary 0.06697091 --viscosity 0.0086466165 --hamaker 6.6428476e-11 --tmax 10 --mfp 0.00046 --usegke  --dropdrop "
if [[ -d $folder ]]; then
echo "$folder already exists, stopping simulation"
exit
fi
mkdir $folder
make
./drop_impact --folder $folder $parameters &> $logfile
