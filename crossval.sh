#!/bin/bash

set -e

if [ -z "$1" ] || [ -z "$2" ];
then
   echo "Usage: [optional params] cv.sh <root name of PLINK file> <model>"
   echo "where model is one of: linear, sqrhinge"
   exit 1
fi

[[ -z "$NUMPROCS" ]] && NUMPROCS=1

D=$(dirname $1)
F=$(basename $1)
ROOT=$D/$F

MODEL=$2

######################################################################
# User modifiable parameters

# Number of cross-validation folds
[[ -z "$NFOLDS" ]] && NFOLDS=10

# Number of cross-validation replications
[[ -z "$NREPS" ]] && NREPS=10

# Number of penalties to look at
[[ -z "$NLAMBDA1" ]] && NLAMBDA1=30

# The smallest L1 penalty, as a proportion of the maximal L1 penalty
# (determined within SparSNP from the data)
[[ -z "$L1MIN" ]] && L1MIN=0.001

# L2 penalty (elastic-net)
[[ -z "$LAMBDA2" ]] && LAMBDA2=0

# Standardise the outputs, only really makes sense for linear regression
# Note that we scale Y globally, not for each cross-validation fold
#[[ -z "$SCALEY" ]] && SCALEY=""

# Fusion penalty
[[ -z "$GAMMA" ]] && GAMMA=0

# By default, return beta on the original scale of the data (before standardising)
UNSCALE=${UNSCALE- "-unscale"}

######################################################################

N=$(cat "$ROOT".fam | wc -l | awk '{print $1, $2}')
P=$(cat "$ROOT".bim | wc -l | awk '{print $1, $2}')
BED=$(realpath "$ROOT".bed)
FAM=$(realpath "$ROOT".fam)
PHENO=$(realpath "$PHENO")
BIM=$(realpath "$ROOT".bim)
SCALE=scale.bin
FOLDIND_CMD="-foldind folds.ind"

# PHENO takes priority over FAM
if ! [ -z "$PHENO" ];
then
   FAM_CMD=""
   PHENO_CMD="-pheno $PHENO"
else
   FAM_CMD="-fam $FAM"
   PHENO_CMD=""
fi



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

echo "N: $N"
echo "P: $N"
echo "NZMAX: $NZMAX"
echo "BED: $BED"
echo "BIM: $BIM"
echo "FAM: $FAM"

DIR=discovery

if ! [ -d "$DIR" ];
then
   mkdir $DIR
fi

cat > $DIR/params.txt<<EOF
ROOT=$(realpath "$ROOT")
FAM=$FAM
PHENO=$PHENO
NFOLDS=$NFOLDS
NREPS=$NREPS
NZMAX=$NZMAX
NLAMBDA1=$NLAMBDA1
LAMBDA2=$LAMBDA2
GAMMA=$GAMMA
MODEL=$MODEL
BETA_SCALED=$UNSCALE
Y_SCALED=$SCALEY
EOF

pushd $DIR
awk '{print $2, $5}' "$BIM" > snps.txt

if [ $NFOLDS -le 1 ];
then
   FOLDIND_CMD=""
fi

[[ -z "$REP_START" ]] && REP_START=1
[[ -z "$REP_END" ]] && REP_END=$NREPS

echo "REP_START: $REP_START"
echo "REP_END: $REP_END"

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
      makefolds -folds folds.txt -ind folds.ind -nfolds $NFOLDS -n $N

      # Get scale of each crossval fold
      scale -bed $BED -n $N -p $P $FOLDIND_CMD
	 
      echo "############# Running training #############"
   
      # Run the model
      $WRAPPER sparsnp -train -model $MODEL -n $N -p $P \
	 -scale $SCALE -bed $BED -nzmax $NZMAX -nl1 $NLAMBDA1 -l1min $L1MIN -v \
	 $FOLDIND_CMD $FAM_CMD $PHENO_CMD -l2 $LAMBDA2 $UNSCALE $SCALEY -gamma $GAMMA
 
      if [ $NFOLDS -gt 1 ]
      then
	 echo "############# Running prediction #############"
	 # Predict for test folds
	 B=$(for((i=0;i<NLAMBDA1;i++)); do printf 'beta.csv.%02d ' $i; done)
	 $WRAPPER sparsnp -predict -model $MODEL -n $N -p $P -v \
	    -bed $BED -betafiles $B \
	    -scale $SCALE \
	    $FOLDIND_CMD $PHENO_CMD
      fi

      popd
   else
      echo "skipping $dir"
   fi
}

export -f run
export NFOLDS N P BED FAM_CMD PHENO_CMD FOLDIND_CMD MODEL
export SCALEY SCALE NZMAX NLAMBDA1 L1MIN LAMBDA2

seq $REP_START $REP_END | xargs -P$NUMPROCS -I{} bash -c "run {}"

popd

