#!/bin/bash

set -e

shopt -s extglob

if [ -z "$1" ];
then
   echo "Usage: validation.sh <root name of PLINK file>"
   exit 1
fi


D=$(dirname $1)
F=$(basename $1)
ROOT=$D/$F

set -u

N=$(cat $ROOT.fam | wc -l)
P=$(cat $ROOT.bim | wc -l)
BIN=$(./realpath "$ROOT".bed)
FAM=$(./realpath "$ROOT".fam)

ORIGWD=$(pwd)
echo "PWD: $PWD"

DIRS=$(ls -d $ORIGWD/discovery/crossval*)

MODEL=$(grep MODEL $ORIGWD/discovery/params.txt | cut -f2 -d'=')

echo "Model for validation: $MODEL"

VDIR=validation
mkdir $VDIR
pushd $VDIR


for dir in $DIRS
do
   pdir=$(basename $dir)
   # don't clobber an existing directory
   if ! [ -d $pdir ];
   then
      mkdir $pdir
      pushd $pdir

      B=$( ls $dir/beta.csv.+([[:digit:]]).+([[:digit:]]) )

      ../../cd -predict -model $MODEL -n $N -p $P -v \
	 -bin $BIN -betafiles $B -fam $FAM -outdir .

      awk '{print $6}' $FAM > y.txt

      popd
   else
      echo "skipping $pdir"
   fi
done

popd

