#!/bin/bash

# Convert SparSNP SNP weights to a PLINK score file
# SNPID refallele score
#
# Requires a recent SparSNP version where the reference allele has been stored

if [ $# -le 2 ];
then
   echo "usage: predict.sh <TARGET ROOT> <non-zero req.> <prevalence>"
   exit 1
fi

set -e

TARGET=$1
NZREQ=$2
PREV=$3

NREPS=10
NFOLDS=10

DIR="discovery"
AVGFILE="avg_weights_opt"

cat > .predict.R <<EOF
nreps <- $NREPS
nfolds <- $NFOLDS
nzreq <- $NZREQ
dir <- "$DIR"
avgfile <- "$AVGFILE"

r <- read.table(sprintf("%s/snps.txt", dir), header=FALSE,
   stringsAsFactors=FALSE)

#lf <- list.files(pattern="^beta.csv.[[:digit:]]+.[[:digit:]]+$",
#   recursive=TRUE, full.names=TRUE)

nz <- lapply(1:nreps, function(i) {
   lapply(1:nfolds, function(j) {
      scan(sprintf("%s/crossval%s/nonzero.csv.%02d", dir, i, j - 1),
	 quiet=TRUE)
   })
})

nzopt <- lapply(nz, sapply, function(n) which.min((n - nzreq)^2))

b <- lapply(1:nreps, function(rep) {
   lapply(1:nfolds, function(fold) {
      w <- nzopt[[rep]][fold]
      f <- sprintf("%s/crossval%s/beta.csv.%02d.%02d",
	 dir, rep, w - 1, fold - 1)
      read.table(f, sep=":", header=FALSE)
   })
})

b <- unlist(b, recursive=FALSE)
b2 <- do.call(rbind, b)
b2m <- with(b2, tapply(V2, V1, mean))

# Average weights at chosen model, over cross-validation reps, 
# including intercept
#write.table(b2m, file=sprintf("%s/%s.txt", dir, avgfile),
#   col.names=FALSE, sep=",", quote=FALSE)
write.table(b2m["0"], file=sprintf("%s/intercept.txt", dir),
   col.names=FALSE, sep=",", quote=FALSE, row.names=FALSE)

# PLINK score file, excluding intercept
b2mp <- b2m[grep("^0$", names(b2m), invert=TRUE)]
snps <- read.table(sprintf("%s/snps.txt", dir))
w <- as.integer(names(b2mp))

d <- data.frame(snps[w,], b2mp)

write.table(d, file=sprintf("%s/%s.score", dir, avgfile),
   quote=FALSE, row.names=FALSE, col.names=FALSE)
EOF

Rscript .predict.R

plink --noweb --bfile $TARGET \
   --score $DIR/$AVGFILE.score \
   --out $(basename $TARGET)

cat > .eval.R <<EOF
library(ROCR)

target <- "$TARGET"
target.base <- "$(basename $TARGET)"
prev <- $PREV
dir <- "$DIR"

fam <- read.table(sprintf("%s.fam", target), header=FALSE, sep="")
prof <- read.table(sprintf("%s.profile", target.base), header=TRUE)
intercept <- scan(sprintf("%s/intercept.txt", dir))

pheno <- fam[, 6]
score <- prof[, 6] + intercept 

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

save(pheno, prof, pred, perf, ppv, npv,
   file=sprintf("%s.RData", target.base))


EOF

Rscript .eval.R

/bin/rm .predict.R .eval.R

