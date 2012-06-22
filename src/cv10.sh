#!/bin/bash

set -e

if [ -z "$1" ] || [ -z "$2" ];
then
   echo "Usage: cv.sh <root name of PLINK file> <model>"
   echo "where model is one of: linear, sqrhinge"
   exit 1
fi

[[ -z "$NUMPROCS" ]] && NUMPROCS=1

D=$(dirname $1)
F=$(basename $1)
ROOT=$D/$F

EXEDIR=~/Code/cd/src/sparsnp

######################################################################
# User modifiable parameters

# number of cross-validation folds
[[ -z "$NFOLDS" ]] && NFOLDS=10

# number of cross-validation replications
[[ -z "$NREPS" ]] && NREPS=10

# number of penalties to look at
[[ -z "$NLAMBDA1" ]] && NLAMBDA1=30

[[ -z "$L1MIN" ]] && L1MIN=0.001

[[ -z "$LAMBDA2" ]] && LAMBDA2=0

MODEL=$2

######################################################################

######################################################################
# Don't change these unless you know what you're doing
N=$(cat "$ROOT".fam | wc -l)
P=$(cat "$ROOT".bim | wc -l)
BED=$($EXEDIR/realpath "$ROOT".bed)
FAM=$($EXEDIR/realpath "$ROOT".fam)
BIM=$($EXEDIR/realpath "$ROOT".bim)
SCALE=scale.bin
FOLDIND="-foldind folds.ind"
######################################################################

# maximum number of non-zero SNPs to consider in model
if [ -z "$NZMAX" ];
then
   if [ $N -lt $P ];
   then
      NZMAX=$N
   else
      NZMAX=$P
   fi
fi

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
ROOT=$ROOT
NFOLDS=$NFOLDS
NREPS=$NREPS
NZMAX=$NZMAX
NLAMBDA1=$NLAMBDA1
LAMBDA2=$LAMBDA2
MODEL=$MODEL
EOF

awk '{print $2, $5}' "$BIM" > snps.txt

if [ $NFOLDS -le 1 ];
then
   FOLDIND=""
fi

[[ -z "$REP_START" ]] && REP_START=1
[[ -z "$REP_END" ]] && REP_END=$NREPS

echo "REP_START: $REP_START"
echo "REP_END: $REP_END"

#for((i=$REP_START;i<=$REP_END;i++))
#do
function run {
   i=$1
   dir="crossval$i"

   echo "dir=$dir"

   # don't clobber an existing directory
   if ! [ -d $dir ];
   then
      mkdir $dir
      pushd $dir

      # Create cross-validation folds
      $EXEDIR/split -folds folds.txt -ind folds.ind -nfolds $NFOLDS -n $N

      # Get scale of each crossval fold
      $EXEDIR/scale -bin $BED -n $N -p $P $FOLDIND 
   
      # Run the model
      $EXEDIR/cd -train -model $MODEL -n $N -p $P \
	 -scale $SCALE -bin $BED -nzmax $NZMAX -nl1 $NLAMBDA1 -l1min $L1MIN -v \
	 $FOLDIND -fam $FAM -l2 $LAMBDA2
 
      # Predict for test folds
      B=$(for((i=0;i<NLAMBDA1;i++)); do printf 'beta.csv.%02d ' $i; done)
      $EXEDIR/cd -predict -model $MODEL -n $N -p $P -v \
	 -bin $BED -betafiles $B \
	 -scale $SCALE \
	 $FOLDIND -fam $FAM

      popd
   else
      echo "skipping $dir"
   fi
#done
}

export -f run
export EXEDIR NFOLDS N P BED FAM FOLDIND MODEL SCALE NZMAX NLAMBDA1 L1MIN LAMBDA2

seq $REP_START $REP_END | xargs -P$NUMPROCS -I{} bash -c "run {}"

popd

