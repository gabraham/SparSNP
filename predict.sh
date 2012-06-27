#!/bin/bash

# Converts SparSNP SNP weights to a PLINK score file
# SNPID refallele score
#
# Requires a recent SparSNP version where the reference allele has been stored

if [ $# -lt 2 ];
then
   echo "usage: predict.sh <TARGET ROOT> <prevalence>"
   exit 1
fi

set -e

TARGET=$1
PREV=$2

DIR="discovery"
OUTDIR="predict"
AVGFILE="avg_weights_opt"

if ! [ -d "$OUTDIR" ];
then
   mkdir "$OUTDIR"
fi

eval $(grep ROOT $DIR/params.txt)

if ! [ -s "$DIR/$AVGFILE.score" ];
then
   echo "Can't find score files, have you run getmodels.R?"
   exit 1
fi

# Get a risk score for the top SNP using logistic regression
plink --noweb --bfile $ROOT \
   --snp $(head -1 $DIR/$AVGFILE.score | cut -f1 -d' ') \
   --logistic --out $DIR/topsnp

awk '{if(NR>1) print $2,$4,log($7)}' \
   $DIR/topsnp.assoc.logistic > $DIR/topsnp.score

# Predict on new data, multi-SNP
for f in $DIR/"$AVGFILE"_path_*.score;
do
   n=$(echo $(basename $f) | sed 's/[avg_weights_opt_path_|\.score]//g')
   plink --noweb --bfile $TARGET \
      --score $f \
      --out $OUTDIR/$(basename $TARGET)_$n
done


