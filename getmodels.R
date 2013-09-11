#!/usr/bin/env Rscript

# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
# 
# Copyright (C) 2011-2012 Gad Abraham and National ICT Australia (NICTA).
# All rights reserved.
# 

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

AVGFILE="avg_weights_opt"

usage <- "usage: getmodels.R nzreq=<num SNPS required> [uni=FALSE] [dir=DIR]"

uni <- FALSE
args <- commandArgs(TRUE)
s <- strsplit(args, "=")
for(m in s) {
   eval(parse(text=sprintf("%s=\"%s\"", m[1], m[2])))
}

if(!exists("nzreq", mode="character")) {
   stop(usage)
}

if(!exists("dir", mode="character")) {
   dir <- "discovery"
}

if(!exists("thresh", mode="character")) {
   thresh <- 0.6
}

nzreq <- as.numeric(nzreq)
if(is.na(nzreq) || as.integer(nzreq) != nzreq || nzreq <= 0) {
   stop(usage)
}

uni <- as.logical(uni)
if(is.na(uni)) {
   stop(usage)
}

library(methods)

# Parse the parameters
param <- scan(sprintf("%s/params.txt", dir), what=character(), quiet=TRUE,
   sep="\n")
s <- sapply(param, strsplit, "=")
for(i in seq(along=s)) {
   s1 <- s[[i]][1]
   if(length(s[[i]]) > 1) {
      s2 <- s[[i]][2]
      m <- as.numeric(s2)
      s2 <- if(is.na(m)) {
         paste("\"", s2, "\"", sep="")
      } else m
   } else {
      s2 <- "\"\""
   }
   x <- paste(s1, "=", s2)
   eval(parse(text=x))
}

if(NTASKS > 1) {
   stop("getmodels currently does not support models with more than 1 task")
}

nreps <- NREPS
nfolds <- NFOLDS
avgfile <- AVGFILE

r <- read.table(sprintf("%s/snps.txt", dir), header=FALSE,
   stringsAsFactors=FALSE)

prefix <- ifelse(uni, "multivar_", "")

nz <- lapply(1:nreps, function(i) {
   lapply(1:nfolds, function(j) {
      scan(sprintf("%s/crossval%s/%snonzero.csv.%02d", dir, i, prefix, j - 1),
	 quiet=TRUE)
   })
})

# Get best model with the closest model size to requested
nzopt <- lapply(nz, sapply, function(n) which.min((n - nzreq)^2))

b <- lapply(1:nreps, function(rep) {
   lapply(1:nfolds, function(fold) {
      w <- nzopt[[rep]][fold]
      f <- sprintf("%s/crossval%s/%sbeta.csv.%02d.%02d",
	 dir, rep, prefix, w - 1, fold - 1)
      r <- try(read.table(f, sep=",", header=FALSE), silent=FALSE)
      if(is(r, "try-error")) {
	 NULL
      } else {
	 r
      }
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
	 pattern=sprintf("^%sbeta.csv.[[:digit:]]+.%02d$", prefix, fold - 1),
	 full.names=TRUE
      )
      lapply(l, function(x) {
	 cat(x, "\n")
	 r <- read.table(x, sep=",", header=FALSE)
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
   w <- names(which(table(r2$V1) / length(r) >= thresh))
   r3 <- r2[r2$V1 %in% w, ]
   with(r3, tapply(V2, V1, function(x) sum(x) / length(r)))
})

# Average weights at chosen model, over cross-validation reps, 
# including intercept
write.table(b2m["0"], file=sprintf("%s/%sintercept.txt", dir, prefix),
   col.names=FALSE, sep=",", quote=FALSE, row.names=FALSE)

# PLINK score file, excluding intercept
b2mp <- b2m[grep("^0$", names(b2m), invert=TRUE)]
b2mp <- b2mp[order(abs(b2mp), decreasing=TRUE)][1:nzreq]
snps <- read.table(sprintf("%s/snps.txt", dir))
w <- as.integer(names(b2mp))

d <- data.frame(snps[w,], b2mp)

write.table(d, file=sprintf("%s/%s%s.score", dir, prefix, avgfile),
   quote=FALSE, row.names=FALSE, col.names=FALSE)

# Write averaged models over the regularisation path, plus intercept for each
# model
for(i in seq(along=bpathm)) {
   x <- bpathm[[i]]
   x <- x[order(abs(x), decreasing=TRUE)]#[1:nzreq]
   w <- as.integer(names(x))
   d <- data.frame(snps[w[w != 0], ], x[w != 0])
   write.table(d, file=sprintf("%s/%s%s_path_%s.score", dir, prefix, avgfile, i),
      quote=FALSE, row.names=FALSE, col.names=FALSE)
   write.table(x["0"], file=sprintf("%s/%sintercept_path_%s.txt", dir, prefix, i),
      col.names=FALSE, sep=",", quote=FALSE, row.names=FALSE)
}


