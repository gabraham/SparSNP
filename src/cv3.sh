#!/bin/bash

set -e


if [ -z "$1" ];
then
   echo "Usage: cv3.sh <root name of PLINK file>"
   exit 1
fi

D=$(dirname $1)
F=$(basename $1)
ROOT=$D/$F


set -u

######################################################################
# User modifiable parameters

# number of cross-validation folds
NFOLDS=3

# number of cross-validation replications
NREPS=3

# maximum number of non-zero SNPs to consider in model
NZMAX=512

# number of penalties to look at
NLAMBDA1=25

# model type, other options: linear
MODEL=sqrhinge

######################################################################

######################################################################
# Don't change these unless you know what you're doing
N=$(cat $ROOT.fam | wc -l)
P=$(cat $ROOT.bim | wc -l)
BIN="$ROOT".bed
PLINK="-plink"
FAM="-fam $ROOT".fam
SCALE=scale.bin
FOLDIND="-foldind folds.ind"
######################################################################

for((i=1;i<=$NREPS;i++))
do
   dir="crossval$i"

   # don't clobber an existing directory
   if ! [ -d $dir ];
   then
      mkdir $dir
      pushd $dir

      # Create cross-validation folds
      ../split -folds folds.txt -ind folds.ind -nfolds $NFOLDS -n $N

      # Get scale of each crossval fold
      ../scale -bin $BIN -n $N -p $P $FOLDIND $PLINK
   
      # Run the model
      ../cd -train -model $MODEL -n $N -p $P \
	 -scale $SCALE -bin $BIN -nzmax $NZMAX -nl1 $NLAMBDA1 -l1min 0.01 -v \
	 $FOLDIND $PLINK $FAM
 
      # Predict for test folds
      B=$(for((i=0;i<=NLAMBDA1;i++)); do printf 'beta.csv.%02d ' $i; done)
      ../cd -predict -model $MODEL -n $N -p $P -v \
	 -bin $BIN -betafiles $B \
	 -scale $SCALE \
	 $FOLDIND $PLINK $FAM

      popd
   else
      echo "skipping $dir"
   fi
done

