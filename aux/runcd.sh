#!/bin/bash

set -u
set -e

DIR_STEM="sim8."
N=30000
P=185805
MODEL=sqrhinge
NZMAX=512
bin=sim.bin
CD="~/Code/cd/src/cd -model $MODEL -bin ../$bin -scale ../scale.bin \
-n $N -p $P -v -maxepochs 10000 -nzmax $NZMAX -nl1 100 -l1min 0.001 \
-thresh 1e-6"
SCALE="~/Code/cd/src/scale -bin $bin -n $N -p $P"
scale="scale.bin"
results="results4"

for((i=10 ; i<=10 ; i++));
do
   dir="$DIR_STEM""$i"
   dir_res="$dir/$results"
   #if ! [ -d "$dir_res" ]; then
    #  mkdir "$dir_res"
      pushd "$dir"
      #if ! [ -f "$scale" ]; then
         eval "$SCALE"   
      #else
      #   echo "skipping scale"
      #fi
      pushd "$results"
      #if ! [ -a "lambda1path.csv" ];
      #then
         eval "time $CD" 2>&1 | tee > log
      #else
      #   echo "skipping cd"
      #fi
      popd
      popd
   #else
   #   echo "skipping $dir_res"
   #fi
done

