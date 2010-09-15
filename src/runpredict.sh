#!/bin/bash

set -u
set -e

shopt -s extglob

DIR_STEM="sim7."
N=1000
P=185805
bin=sim.bin
enc="-encoded"
scale="scale.bin"
results="results"

for((i=1 ; i<=10 ; i++));
do
   dir="$DIR_STEM""$i"
   while ! [[ -d "$dir" && -a "$dir/$bin" ]];
   do
      echo "$dir/$bin doesn't exist yet, sleeping..."
      sleep 300
   done
   dir_res="$dir/$results"
   if ! [ -d "$dir_res" ]; then
      continue
   fi
   pushd "$dir/$results"
   
   # predict the phenotype
   CD="~/Code/cd/src/cd -predict -model sqrhinge \
-f ../$bin -scale ../scale.bin \
-n $N -p $P -v -betafiles $(ls -v beta.csv.+([0-9])) \
$enc"
   eval "$CD"

   # unscale the coefs
   for((k=0 ; k <= 30 ; k++));
   do
      UNSC="~/Code/cd/src/scale -unscale -p $P -scale ../scale.bin \
-betafile beta.csv.$k -enc"
      eval "$UNSC"
   done

   popd
done

