#!/bin/bash

# This script MUST be called after crossval.sh has been called

set -e


if [ -z "$1" ] ; 
then
   echo "Usage: crossvaluni.sh <root name of PLINK file>"
   exit 1
fi

[[ -z "$NUMPROCS" ]] && NUMPROCS=1

D=$(dirname $1)
F=$(basename $1)
ROOT=$D/$F

EXEDIR=~/Code/cd/src/sparsnp
MODEL="logistic"

# Don't change these unless you know what you're doing
N=$(cat "$ROOT".fam | wc -l)
P=$(cat "$ROOT".bim | wc -l)
BED=$($EXEDIR/realpath "$ROOT".bed)
FAM=$($EXEDIR/realpath "$ROOT".fam)
BIM=$($EXEDIR/realpath "$ROOT".bim)

echo "BED: $BED"
echo "BIM: $BIM"
echo "FAM: $FAM"

DIR=discovery

if ! [ -d "$DIR" ];
then
   echo "Dir $DIR doesn't exist, did you run crossval.sh first?"
   exit 1
fi

pushd $DIR

eval $(grep NREPS params.txt)

[[ -z "$NZMAX" ]] && NZMAX=20000
[[ -z "$REP_START" ]] && REP_START=1
[[ -z "$REP_END" ]] && REP_END=$NREPS

function run {
   i=$1
   dir="crossval$i"

   echo "dir=$dir"

   # don't clobber an existing directory
   if ! [ -d $dir ];
   then
      echo "skipping $dir, doesn't exist"
   else
      pushd $dir

      $EXEDIR/univariable -train -model $MODEL \
	 -bin $BED -fam $FAM -n $N -p $P -nzmax $NZMAX \
         -foldind folds.ind -v

      nthresh=$(cat multivar_nonzero.csv.00 | wc -l)
      tfiles=$(for((i=0;i<$nthresh;i++)); \
            do printf "multivar_beta.csv.%02d " $i; \
            done)

      $EXEDIR/univariable -predict -model logistic \
	 -bin $BED -fam $FAM -n $N -p $P \
         -foldind folds.ind -v -predict -betafiles $tfiles

      popd
   fi
}

export -f run
export EXEDIR NFOLDS N P BED FAM MODEL NZMAX 

seq $REP_START $REP_END | xargs -P$NUMPROCS -I{} bash -c "run {}"

popd

