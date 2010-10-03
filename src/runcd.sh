#!/bin/bash

set -u
set -e

DIR_STEM="sim7."
N=10000
P=185805
bin=sim.bin
CD="~/Code/cd/src/cd -model sqrhinge -bin ../$bin -scale ../scale.bin \
-n $N -p $P -v -maxepochs 10000 -nzmax 250 -nl1 100"
SCALE="~/Code/cd/src/scale -bin $bin -n $N -p $P"
scale="scale.bin"
results="results4"

for((i=21 ; i<=21 ; i++));
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

