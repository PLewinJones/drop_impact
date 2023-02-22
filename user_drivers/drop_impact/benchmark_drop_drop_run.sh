#!/bin/bash 

folder=RESLT_benchmark_drop_drop 
logfile=$folder"/oomph-out.txt"
parameters="--reynolds 29.767272 --gravity 0 --capilliary 0.039829736 --viscosity 0.0085855263 --hamaker 8.9467964e-11 --tmax 10 --mfp 0.00041169451 --usegke  --dropdrop "
if [[ -d $folder ]]; then
echo "$folder already exists, stopping simulation"
exit
fi
mkdir $folder
make
./drop_impact --folder $folder $parameters &> $logfile
