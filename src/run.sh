#!/bin/bash

N=1000
P=185805
NFOLDS=10
NZMAX=500

~/Code/cd/src/split -folds folds.txt -ind folds.ind -nfolds $NFOLDS -n $N

~/Code/cd/src/scale -bin ../sim.bin -n $N -p $P -foldind folds.ind

time ~/Code/cd/src/cd -train -model sqrhinge -n $N -p $P -v \
-foldind folds.ind -scale scale.bin -bin ../sim.bin -nzmax $NZMAX \
-maxepochs 5000

B=$(for((i=0;i<=50;i++)); do printf 'beta.csv.%02d ' $i; done)
~/Code/cd/src/cd -predict -model sqrhinge -n $N -p $P -v \
-foldind folds.ind -bin ../sim.bin -betafiles $B -scale scale.bin


