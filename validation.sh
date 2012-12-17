#!/bin/bash

set -e

shopt -s extglob

if [ -z "$1" ];
then
   echo "Usage: validation.sh <root name of PLINK file>"
   exit 1
fi

[[ -z "$NUMPROCS" ]] && NUMPROCS=1

D=$(dirname $1)
F=$(basename $1)
ROOT=$D/$F

set -u

N=$(cat $ROOT.fam | wc -l)
P=$(cat $ROOT.bim | wc -l)
BED=$(realpath "$ROOT".bed)
FAM=$(realpath "$ROOT".fam)

ORIGWD=$(pwd)
echo "PWD: $PWD"

DIRS=$(ls -d $ORIGWD/discovery/crossval*)

MODEL=$(grep MODEL $ORIGWD/discovery/params.txt | cut -f2 -d'=')
NFOLDS=$(grep NFOLDS $ORIGWD/discovery/params.txt | cut -f2 -d'=')

echo "Model for validation: $MODEL"

VDIR=validation
mkdir $VDIR
pushd $VDIR


function run {
   dir=$1
   pdir=$(basename $dir)
   # don't clobber an existing directory
   if ! [ -d $pdir ];
   then
      mkdir $pdir
      pushd $pdir

      shopt -s extglob

      B=` ls $dir/beta.csv.+([[:digit:]]).+([[:digit:]]) `

      echo "Found files $B"

      sparsnp -predict -model $MODEL -n $N -p $P -v \
         -bed $BED -betafiles $B -fam $FAM -outdir .

      awk '{print $6}' $FAM > y.txt

      if [ -a "$dir/multivar_nonzero.csv.00" ];
      then
         nthresh=$(cat $dir/multivar_nonzero.csv.00 | wc -l)
   
         if [ $nthresh -gt 0 ];
         then
   	 #tfiles=$(for((i=0;i<$nthresh;i++)); \
            #   do printf "$dir/multivar_beta.csv.%02d " $i; \
            #   done)
   	 B=$( ls $dir/multivar_beta.csv.+([[:digit:]]).+([[:digit:]]) )
   
   	 # Predict on entire validation dataset using each model
   	 univariable -predict -model logistic \
   	    -bin $BED -fam $FAM -n $N -p $P \
   	    -betafiles $B \
   	    -outdir .
   
         fi
      fi

      popd
   else
      echo "skipping $pdir"
   fi
}

export -f run
export MODEL N P BED B FAM

echo $DIRS | tr ' ' '\n' | xargs -P$NUMPROCS -I{} bash -c "run {}"

popd

