#!/bin/bash

# Converts SparSNP SNP weights to a PLINK score file
# SNPID refallele score
#
# Requires a recent SparSNP version where the reference allele has been stored

if [ $# -lt 2 ];
then
   #echo "usage: predict.sh <TARGET ROOT> <non-zero req.> <prevalence>"
   echo "usage: predict.sh <TARGET ROOT> <prevalence>"
   exit 1
fi

set -e

TARGET=$1
#NZREQ=$2
#PREV=$3
PREV=$2

DIR="discovery"
OUTDIR="predict"
AVGFILE="avg_weights_opt"

if ! [ -d "$OUTDIR" ];
then
   mkdir "$OUTDIR"
fi

#eval $(grep NFOLDS $DIR/params.txt)
#eval $(grep NREPS $DIR/params.txt)
eval $(grep ROOT $DIR/params.txt)

## Averages models along the regularisation path, across cross-validation
## replications, and outputs these models plus the best model,
## in the form of a PLINK score file
#cat > .predict.R <<EOF
#nreps <- $NREPS
#nfolds <- $NFOLDS
#nzreq <- $NZREQ
#dir <- "$DIR"
#avgfile <- "$AVGFILE"
#
#library(Hmisc)
#
#r <- read.table(sprintf("%s/snps.txt", dir), header=FALSE,
#   stringsAsFactors=FALSE)
#
#nz <- lapply(1:nreps, function(i) {
#   lapply(1:nfolds, function(j) {
#      scan(sprintf("%s/crossval%s/nonzero.csv.%02d", dir, i, j - 1),
#	 quiet=TRUE)
#   })
#})
#
#nzopt <- lapply(nz, sapply, function(n) which.min((n - nzreq)^2))
#
#b <- lapply(1:nreps, function(rep) {
#   lapply(1:nfolds, function(fold) {
#      w <- nzopt[[rep]][fold]
#      f <- sprintf("%s/crossval%s/beta.csv.%02d.%02d",
#	 dir, rep, w - 1, fold - 1)
#      read.table(f, sep=":", header=FALSE)
#   })
#})
#
#b <- unlist(b, recursive=FALSE)
#b2 <- do.call(rbind, b)
#b2m <- with(b2, tapply(V2, V1, function(x) sum(x) / length(b)))
#
## Evaluate the entire path of solutions, averaging over all folds and
## replications and clusters of solutions with roughly same size.
## We don't check cluster size, just take the kth cluster across all
## replications, as with a fine enough penalty grid model sizes will be
## similar.
#bpath <- lapply(1:nreps, function(rep) {
#   lapply(1:nfolds, function(fold) {
#      l <- list.files(
#	 path=sprintf("%s/crossval%s", dir, rep),
#	 pattern=sprintf("^beta.csv.[[:digit:]]+.%02d$", fold - 1),
#	 full.names=TRUE
#      )
#      lapply(l, function(x) {
#	 r <- read.table(x, sep=":", header=FALSE)
#	 r\$NonZero <- nrow(r)
#	 r
#      })
#   })
#})
#bpath2 <- unlist(bpath, recursive=FALSE)
#m <- min(sapply(bpath2, length))
#bpath3 <- lapply(bpath2, head, n=m)
#bpathm <- lapply(1:m, function(i) {
#   r <- lapply(bpath3, function(x) x[[i]])
#   r2 <- do.call(rbind, r) 
#   with(r2, tapply(V2, V1, function(x) sum(x) / length(r)))
#})
#
## Average weights at chosen model, over cross-validation reps, 
## including intercept
#write.table(b2m["0"], file=sprintf("%s/intercept.txt", dir),
#   col.names=FALSE, sep=",", quote=FALSE, row.names=FALSE)
#
## PLINK score file, excluding intercept
#b2mp <- b2m[grep("^0$", names(b2m), invert=TRUE)]
#b2mp <- b2mp[order(abs(b2mp), decreasing=TRUE)][1:nzreq]
#snps <- read.table(sprintf("%s/snps.txt", dir))
#w <- as.integer(names(b2mp))
#
#d <- data.frame(snps[w,], b2mp)
#
#write.table(d, file=sprintf("%s/%s.score", dir, avgfile),
#   quote=FALSE, row.names=FALSE, col.names=FALSE)
#
## Write averaged models over the regularisation path
#for(i in seq(along=bpathm))
#{
#   x <- bpathm[[i]]
#   x <- x[order(abs(x), decreasing=TRUE)]
#   w <- as.integer(names(x))
#   d <- data.frame(snps[w[w != 0], ], x[w != 0])
#   write.table(d, file=sprintf("%s/%s_path_%s.score", dir, avgfile, i),
#      quote=FALSE, row.names=FALSE, col.names=FALSE)
#   write.table(x["0"], file=sprintf("%s/intercept_path_%s.txt", dir, i),
#      col.names=FALSE, sep=",", quote=FALSE, row.names=FALSE)
#}
#
#EOF
#
#Rscript .predict.R

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

# Predict on new data, single SNP
plink --noweb --bfile $TARGET \
   --score $DIR/topsnp.score \
   --out $OUTDIR/$(basename $TARGET)_top1

cat > .eval.R <<EOF
library(ROCR)

r <- commandArgs(TRUE)
target <- r[1]
target.base <- r[2]
famf <- r[3]
prev <- as.numeric(r[4])
dir <- r[5]

fam <- read.table(famf, header=FALSE, sep="")
prof <- read.table(sprintf("%s.profile", target.base), header=TRUE)
intercept <- if(!grepl("top1", target)) {
   scan(sprintf("%s/intercept.txt", dir))
} else 0

pheno <- fam[, 6]
score <- prof[, 6] + intercept 

# Prevent ties in data that results in fewer cutoffs than data points
dither <- function(x, mind=1e-20, maxd=1e-8)
{
   x + runif(length(x), min=mind, max=maxd)
}
score <- dither(score)

pred <- prediction(labels=fam[,6] - 1, predictions=score)
perf <- performance(pred, "sens", "spec")
sens <- perf@y.values
spec <- perf@x.values
cutoffs <- pred@cutoffs
  
ppv <- sapply(1:length(cutoffs), function(j) {
   (sens[[j]] * prev) / (
      sens[[j]] * prev + (1 - spec[[j]]) * (1 - prev))
})
npv <- sapply(1:length(cutoffs), function(j) {
   (spec[[j]] * (1 - prev)) / (
      (1 - sens[[j]]) * prev + spec[[j]] * (1 - prev))
})

save(pheno, score, intercept, prev, prof, pred, perf, sens, spec,
   ppv, npv, file=sprintf("%s.RData", target.base))

EOF

FAMF=$TARGET.fam
Rscript .eval.R $TARGET $(basename $TARGET) $FAMF $PREV $DIR
Rscript .eval.R "$TARGET"_top1 $(basename "$TARGET"_top1) $FAMF $PREV $DIR

#/bin/rm .predict.R .eval.R

