#!/usr/bin/env Rscript

# getmodels.R:
# Averages models along the regularisation path, across 
# cross-validation replications, and outputs these models
# plus the best model, in the form of a PLINK score file
# 
#
# Note: we assume that the solution path produces roughly the same number
# of SNPs in the model across all cross-validation replications. When sample
# sizes are small, or the signal is very weak, this may not hold.
#
# Also, if you ask for a model size that is far from the model sizes that have
# been run, then you won't get back a prediction based on that model size,
# instead you'll get a prediction from a the closest model that could be
# found.

options(warn=-1)

DIR="discovery"
AVGFILE="avg_weights_opt"

usage <- "usage: getmodels.R nzreq=<num SNPS required>"

args <- commandArgs(TRUE)
s <- strsplit(args, "=")
for(m in s) {
   eval(parse(text=sprintf("%s=\"%s\"", m[1], m[2])))
}

if(!exists("nzreq", mode="character")) {
   stop(usage)
}

nzreq <- as.numeric(nzreq)
if(is.na(nzreq) || as.integer(nzreq) != nzreq || nzreq <= 0) {
   stop(usage)
}

library(methods)

# Parse the parameters
param <- scan("discovery/params.txt", what=character(), quiet=TRUE)
s <- sapply(param, strsplit, "=")
for(i in seq(along=s)) {
   s1 <- s[[i]][1]
   s2 <- s[[i]][2]
   m <- as.numeric(s2)
   s2 <- if(is.na(m)) {
      paste("\"", s2, "\"", sep="")
   } else m
   x <- paste(s1, "=", s2)
   eval(parse(text=x))
}

nreps <- NREPS
nfolds <- NFOLDS
dir <- DIR
avgfile <- AVGFILE

r <- read.table(sprintf("%s/snps.txt", dir), header=FALSE,
   stringsAsFactors=FALSE)

nz <- lapply(1:nreps, function(i) {
   lapply(1:nfolds, function(j) {
      scan(sprintf("%s/crossval%s/nonzero.csv.%02d", dir, i, j - 1),
	 quiet=TRUE)
   })
})

# Get best model with the closest model size to requested
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
b2m <- with(b2, tapply(V2, V1, function(x) sum(x) / length(b)))

# Evaluate the entire path of solutions, averaging over all folds and
# replications and clusters of solutions with roughly same size.
# We don't check cluster size, just take the kth cluster across all
# replications, as with a fine enough penalty grid model sizes will be
# similar.
bpath <- lapply(1:nreps, function(rep) {
   lapply(1:nfolds, function(fold) {
      l <- list.files(
	 path=sprintf("%s/crossval%s", dir, rep),
	 pattern=sprintf("^beta.csv.[[:digit:]]+.%02d$", fold - 1),
	 full.names=TRUE
      )
      lapply(l, function(x) {
	 r <- read.table(x, sep=":", header=FALSE)
	 r$NonZero <- nrow(r)
	 r
      })
   })
})

bpath2 <- unlist(bpath, recursive=FALSE)
m <- min(sapply(bpath2, length))
bpath3 <- lapply(bpath2, head, n=m)
bpathm <- lapply(1:m, function(i) {
   r <- lapply(bpath3, function(x) x[[i]])
   r2 <- do.call(rbind, r) 
   with(r2, tapply(V2, V1, function(x) sum(x) / length(r)))
})

# Average weights at chosen model, over cross-validation reps, 
# including intercept
write.table(b2m["0"], file=sprintf("%s/intercept.txt", dir),
   col.names=FALSE, sep=",", quote=FALSE, row.names=FALSE)

# PLINK score file, excluding intercept
b2mp <- b2m[grep("^0$", names(b2m), invert=TRUE)]
b2mp <- b2mp[order(abs(b2mp), decreasing=TRUE)][1:nzreq]
snps <- read.table(sprintf("%s/snps.txt", dir))
w <- as.integer(names(b2mp))

d <- data.frame(snps[w,], b2mp)

write.table(d, file=sprintf("%s/%s.score", dir, avgfile),
   quote=FALSE, row.names=FALSE, col.names=FALSE)

# Write averaged models over the regularisation path, plus intercept for each
# model
for(i in seq(along=bpathm)) {
   x <- bpathm[[i]]
   x <- x[order(abs(x), decreasing=TRUE)][1:nzreq]
   w <- as.integer(names(x))
   d <- data.frame(snps[w[w != 0], ], x[w != 0])
   write.table(d, file=sprintf("%s/%s_path_%s.score", dir, avgfile, i),
      quote=FALSE, row.names=FALSE, col.names=FALSE)
   write.table(x["0"], file=sprintf("%s/intercept_path_%s.txt", dir, i),
      col.names=FALSE, sep=",", quote=FALSE, row.names=FALSE)
}


