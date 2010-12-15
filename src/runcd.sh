#!/bin/bash

set -u
set -e

DIR_STEM="sim4."
N=200
P=100
MODEL=linear
NZMAX=10
bin=sim.bin
CD="~/Code/cd/src/cd -model $MODEL -bin ../$bin -scale ../scale.bin \
-n $N -p $P -v -maxepochs 10000 -nzmax $NZMAX -nl1 100 -l1min 0.001 \
-thresh 1e-6"
SCALE="~/Code/cd/src/scale -bin $bin -n $N -p $P"
scale="scale.bin"
results="results4"

for((i=3 ; i<=3 ; i++));
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
   else
      echo "skipping scale"
   fi
   pushd "$results"
   if ! [ -a "lambda1path.csv" ];
   then
      eval "time $CD" 2>&1 | tee > log
   else
      echo "skipping cd"
   fi
   popd
   popd
done

