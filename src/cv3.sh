#!/bin/bash

set -e


if [ -z "$1" ] || [ -z "$2" ];
then
   echo "Usage: cv3.sh <root name of PLINK file> <model>"
   echo "where model is one of: linear, sqrhinge"
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
NREPS=10

# maximum number of non-zero SNPs to consider in model
NZMAX=1024
#NZMAX=50

# number of penalties to look at
NLAMBDA1=20
#NLAMBDA1=30

L1MIN=0.01

MODEL=$2

######################################################################

######################################################################
# Don't change these unless you know what you're doing
N=$(cat "$ROOT".fam | wc -l)
P=$(cat "$ROOT".bim | wc -l)
BED=$(./realpath "$ROOT".bed)
FAM=$(./realpath "$ROOT".fam)
BIM=$(./realpath "$ROOT".bim)
SCALE=scale.bin
FOLDIND="-foldind folds.ind"
######################################################################

echo "BED: $BED"
echo "BIM: $BIM"
echo "FAM: $FAM"

DIR=discovery

if ! [ -d "$DIR" ];
then
   mkdir $DIR
fi
pushd $DIR

cat >params.txt<<EOF
NFOLDS=$NFOLDS
NREPS=$NREPS
NZMAX=$NZMAX
NLAMBDA1=$NLAMBDA1
MODEL=$MODEL
EOF

awk '{print $2}' "$BIM" > snps.txt


for((i=1;i<=$NREPS;i++))
do
   dir="crossval$i"

   # don't clobber an existing directory
   if ! [ -d $dir ];
   then
      mkdir $dir
      pushd $dir

      # Create cross-validation folds
      ../../split -folds folds.txt -ind folds.ind -nfolds $NFOLDS -n $N

      # Get scale of each crossval fold
      ../../scale -bin $BED -n $N -p $P $FOLDIND 
   
      # Run the model
      ../../cd -train -model $MODEL -n $N -p $P \
	 -scale $SCALE -bin $BED -nzmax $NZMAX -nl1 $NLAMBDA1 -l1min $L1MIN -v \
	 $FOLDIND -fam $FAM
 
      # Predict for test folds
      B=$(for((i=0;i<NLAMBDA1;i++)); do printf 'beta.csv.%02d ' $i; done)
      ../../cd -predict -model $MODEL -n $N -p $P -v \
	 -bin $BED -betafiles $B \
	 -scale $SCALE \
	 $FOLDIND -fam $FAM

      popd
   else
      echo "skipping $dir"
   fi
done

popd

