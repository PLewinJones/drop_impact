#!/bin/bash 

folder=RESLT_benchmark_drop_wall 
logfile=$folder"/oomph-out.txt"
parameters="--reynolds 32.75844 --gravity 3.2801569 --capilliary 0.065989848 --viscosity 0.003654 --hamaker 3.1714588e-12 --tmax 10 --mfp 0.0001 --usegke "
if [[ -d $folder ]]; then
echo "$folder already exists, stopping simulation"
exit
fi
mkdir $folder
make
./drop_impact --folder $folder $parameters &> $logfile
