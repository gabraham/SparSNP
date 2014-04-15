#!/bin/bash

#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
# 
# Copyright (C) 2011-2012 Gad Abraham and National ICT Australia (NICTA).
# All rights reserved.
#

set -e

if [ -z "$1" ] || [ -z "$2" ];
then
   echo "Usage: [optional params] crossval.sh <root name of PLINK file> <model>"
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

# L2 fusion penalty (FMPR)
[[ -z "$GAMMA" ]] && GAMMA=0

# Standardise the outputs, only really makes sense for linear regression
# Note that we scale Y globally, not for each cross-validation fold
! [[ -z "$SCALEY" ]] && SCALEY="-scaley"

# By default, return beta on the original scale of the data (before standardising)
UNSCALE=${UNSCALE-"-unscale"}

[[ -z "$VERBOSE" ]] && VERBOSE="-v"

[[ -z "$DIR" ]] && DIR="discovery"

! [[ -z "$LAMBDA1PATH" ]] && LAMBDA1PATH="-lambda1pathfile_input $(realpath $LAMBDA1PATH)"

######################################################################

N=$(cat "$ROOT".fam | wc -l | awk '{print $1, $2}')
P=$(cat "$ROOT".bim | wc -l | awk '{print $1, $2}')
BED=$(realpath "$ROOT".bed)
FAM=$(realpath "$ROOT".fam)
! [[ -z "$PHENO" ]] && PHENO=$(realpath "$PHENO")
BIM=$(realpath "$ROOT".bim)
SCALE=scale.bin
FOLDIND_CMD="-foldind folds.ind"

# PHENO takes priority over FAM
if ! [ -z "$PHENO" ];
then
   FAM_CMD=""
   PHENO_CMD="-pheno $PHENO"
   NCOL=$(head -1 $PHENO | wc -w)
   NTASKS=$((NCOL-2))
else
   FAM_CMD="-fam $FAM"
   PHENO_CMD=""
   NTASKS=1
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

if ! [ -d "$DIR" ];
then
   mkdir $DIR
fi

cat > $DIR/params.txt<<EOF
ROOT=$(echo $BED | sed 's/\.bed$//g')
BED=$BED
SCALE=$SCALE
N=$N
P=$P
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
NTASKS=$NTASKS
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
   if ! [ -d $dir ] || [ "$CLOBBER" ];
   then
      mkdir $dir
      pushd $dir

      # Get scale of each crossval fold
      if [ "$USEFOLDS" ]
      then
	 echo "Using old folds and scale files."
      else
	 # Create cross-validation folds
	 makefolds -folds folds.txt -ind folds.ind -nfolds $NFOLDS -n $N
	 scale -bed $BED -n $N -p $P $FOLDIND_CMD
      fi
	 
      echo "############# Running training #############"
   
      # Run the model
      $WRAPPER sparsnp -train -model $MODEL -n $N -p $P \
	 -scale $SCALE -bed $BED -nzmax $NZMAX -nl1 $NLAMBDA1 -l1min $L1MIN \
	 $FOLDIND_CMD $FAM_CMD $PHENO_CMD -l2 $LAMBDA2 \
	 -gamma $GAMMA $UNSCALE $SCALEY $VERBOSE $LAMBDA1PATH
 
      if [ $NFOLDS -gt 1 ]
      then
	 echo "############# Running prediction #############"
	 # Predict for test folds
	 B=$(for((i=0;i<NLAMBDA1;i++)); do printf 'beta.csv.%02d ' $i; done)
	 $WRAPPER sparsnp -predict -model $MODEL -n $N -p $P \
	    -bed $BED -betafiles $B \
	    -scale $SCALE $SCALEY \
	    $FOLDIND_CMD $FAM_CMD $PHENO_CMD \
	    $VERBOSE
      fi

      popd
   else
      echo "skipping $dir"
   fi
}

export -f run
export NFOLDS N P BED FAM_CMD PHENO_CMD FOLDIND_CMD MODEL
export UNSCALE SCALEY SCALE NZMAX NLAMBDA1 L1MIN LAMBDA2 GAMMA
export VERBOSE

seq $REP_START $REP_END | xargs -P$NUMPROCS -I{} bash -c "run {}"

popd

