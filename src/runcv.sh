#!/bin/bash

N=500
#P=185805
P=1000
NFOLDS=10
NZMAX=1000
MODEL=sqrhinge

~/Code/cd/src/unpack -bin ../sim.bin -out ../sim -n $N -p $P

~/Code/cd/src/split -folds folds.txt -ind folds.ind -nfolds $NFOLDS -n $N

~/Code/cd/src/scale -bin ../sim.bin -n $N -p $P -foldind folds.ind

time ~/Code/cd/src/cd -train -model $MODEL -n $N -p $P \
-scale scale.bin -bin ../sim.bin -nzmax $NZMAX \
-maxepochs 5000 \
-foldind folds.ind

B=$(for((i=0;i<=100;i++)); do printf 'beta.csv.%02d ' $i; done)
~/Code/cd/src/cd -predict -model $MODEL -n $N -p $P -v \
-bin ../sim.bin -betafiles $B -scale scale.bin \
-foldind folds.ind


