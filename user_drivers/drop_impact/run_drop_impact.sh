#!/bin/bash 

folder=RESLT_example 
logfile=$folder"/oomph-out.txt"
parameters="--reynolds 26.856203 --gravity 0 --capilliary 0.040150943 --viscosity 0.0086466165 --hamaker 1.1080127e-10 --tmax 10 --mfp 0.00046 --usegke  --dropdrop "
if [[ -d $folder ]]; then
echo "$folder already exists, stopping simulation"
exit
fi
mkdir $folder
make
./drop_impact --folder $folder $parameters &> $logfile
