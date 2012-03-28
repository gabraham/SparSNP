#!/bin/bash

set -e


if [ -z "$1" ] ; 
then
   echo "Usage: cv3univar.sh <root name of PLINK file>"
   exit 1
fi

D=$(dirname $1)
F=$(basename $1)
ROOT=$D/$F


set -u

MODEL="logistic"

# Don't change these unless you know what you're doing
N=$(cat "$ROOT".fam | wc -l)
P=$(cat "$ROOT".bim | wc -l)
BED=$(./realpath "$ROOT".bed)
FAM=$(./realpath "$ROOT".fam)
BIM=$(./realpath "$ROOT".bim)

echo "BED: $BED"
echo "BIM: $BIM"
echo "FAM: $FAM"

DIR=discovery

if ! [ -d "$DIR" ];
then
   echo "Dir $DIR doesn't exist, did you run cv3.sh first?"
   exit 1
fi
pushd $DIR

eval $(grep NREPS params.txt)

for((i=1;i<=$NREPS;i++))
do
   dir="crossval$i"

   # don't clobber an existing directory
   if ! [ -d $dir ];
   then
      echo "skipping $dir, doesn't exist"
   else
      pushd $dir

      ../../univariable -train -model logistic \
	 -bin $BED -fam $FAM -n $N -p $P \
         -foldind folds.ind -v

      nthresh=$(cat multivar_nonzero.csv.00 | wc -l)
      tfiles=$(for((i=0;i<$nthresh;i++)); \
            do printf "multivar_beta.csv.%02d " $i; \
            done)

      ../../univariable -predict -model logistic \
	 -bin $BED -fam $FAM -n $N -p $P \
         -foldind folds.ind -v -predict -betafiles $tfiles

      popd
   fi
done

