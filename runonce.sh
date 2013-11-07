#!/bin/bash

set -e

if [[ $# -lt 3 ]];
then
   echo "Usage: runonce.sh <root name of PLINK file> <model> <lambda1>"
   echo "where model is one of: linear, sqrhinge and lambda1 is a positive number"
   exit 1
fi

ROOT=$1
MODEL=$2
L1=$3

N=$(cat "$ROOT".fam | wc -l | awk '{print $1, $2}')
P=$(cat "$ROOT".bim | wc -l | awk '{print $1, $2}')

echo $N $P

makefolds -bed ${ROOT}.bed -n $N -p $P \
   -nfolds 1 -ind folds.ind -folds folds.txt 
scale -bed ${ROOT}.bed -n $N -p $P -foldind folds.ind
sparsnp -train -model $MODEL -n $N -p $P -foldind folds.ind -l1 $L1 \
   -v -bed ${ROOT}.bed -fam ${ROOT}.fam -scale scale.bin.00

