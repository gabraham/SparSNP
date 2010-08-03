#!/bin/bash

set -u
set -e

DIR_STEM="sim6."
N=2000
P=185805
CD="~/Code/cd/src/cd -model sqrhinge -f ../sim.bin.t -scale ../scale.bin \
-n $N -p $P -v -beta beta_sqrhinge.csv"
SCALE="~/Code/cd/src/scale -fin sim.bin.t -n $N -p $P"
scale="scale.bin"
results="results"

for i in $(seq 1 10);
do
   dir="$DIR_STEM""$i"
   dir_res="$dir/$results"
   if ! [ -d "$dir_res" ]; then
      mkdir "$dir_res"
   fi
   pushd "$dir"
   if ! [ -f "$scale" ]; then
      eval "$SCALE"   
   fi
   pushd "$results"
   eval "time $CD" 2>&1 | tee > log
   popd
   popd
done

