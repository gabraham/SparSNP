#!/bin/bash

set -u
set -e

DIR_STEM="sim8."
N=30000
P=185805
bin=sim.bin
enc="-encoded"
CD="~/Code/cd/src/cd -model sqrhinge -f ../$bin -scale ../scale.bin \
-n $N -p $P -v -betafiles beta_sqrhinge.csv -maxepochs 1000 -nzmax 1000 \
-nl1 50 $enc"
SCALE="~/Code/cd/src/scale -fin $bin -n $N -p $P $enc"
scale="scale.bin"
results="results"

for i in $(seq 1 2);
do
   dir="$DIR_STEM""$i"
   while ! [[ -d "$dir" && -a "$dir/$bin" ]];
   do
      echo "$dir/$bin doesn't exist yet, sleeping..."
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

