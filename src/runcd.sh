#!/bin/bash

set -u
set -e

DIR_STEM="sim7."
N=10000
P=185805
CD="~/Code/cd/src/cd -model sqrhinge -f ../sim.bin.t -scale ../scale.bin \
-n $N -p $P -v -betafiles beta_sqrhinge.csv -epochs 1000 -nzmax 1000"
SCALE="~/Code/cd/src/scale -fin sim.bin.t -n $N -p $P"
scale="scale.bin"
results="results"

for i in $(seq 21 30);
do
   dir="$DIR_STEM""$i"
   while ! [[ -d "$dir" && -a "$dir/sim.bin.t" ]];
   do
      echo "Dir $dir doesn't exist yet, sleeping..."
      sleep 300
   done
   dir_res="$dir/$results"
   if ! [ -d "$dir_res" ]; then
      mkdir "$dir_res"
   fi
   pushd "$dir"
   if ! [ -f "$scale" ]; then
      eval "$SCALE"   
   fi
   pushd "$results"
   if ! [ -a "lambda1path.csv" ];
   then
      eval "time $CD" 2>&1 | tee > log
   fi
   popd
   popd
done

