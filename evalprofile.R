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

usage <- paste("usage: evalprofile.R",
   "model=<MODEL>",
   "[prev=<PREVALENCE>]",
   "[name=results]",
   "[indir=discovery]",
   "[outdir=predict]"
)

args <- commandArgs(TRUE)
s <- strsplit(args, "=")
for(m in s) {
   eval(parse(text=sprintf("%s=\"%s\"", m[1], m[2])))
}

if(!exists("outdir", mode="character")) {
   outdir <- "predict"
}

if(!exists("indir", mode="character")) {
   indir <- "discovery"
}

if(!exists("name")) {
   name <- "results"
}

if(!exists("prev")) {
   prev <- NULL
   cat("warning: prevalence not supplied, PPV/NPV will not be calculated\n")
} else {
   prev <- as.numeric(prev)
   if(is.na(prev) || prev < 0 || prev > 1) {
      stop("invalid prevalence: ", prev)
   }
}

library(ROCR)

cat("indir:", indir, "\n")
cat("outdir:", outdir, "\n")

# Prevent ties in data that results in fewer cutoffs than data points
dither <- function(x, mind=1e-20, maxd=1e-8)
{
   x + runif(length(x), min=mind, max=maxd)
}

lf <- list.files(path=outdir, pattern="profile$", full.names=TRUE)
cat("found", length(lf), "profile files\n")
nums <- sapply(sapply(strsplit(lf, "\\.profile"), strsplit, split="_"), tail, n=1)

res <- lapply(seq(along=lf), function(i) {
   cat("reading", lf[i], "\n")
   prof <- read.table(lf[i], header=TRUE)
   intf <- sprintf("%s/intercept_path_%s.txt", indir, nums[i]) 
   cat("reading intercept file", intf, "\n")
   intercept <- scan(intf, quiet=TRUE)
   
   # PLINK divides the predicted score by the number of SNPs, we don't want
   # that to we multiply to get original score
   score <- dither(prof$SCORE * prof$CNT + intercept)
   pred <- prediction(labels=prof$PHENO, predictions=score)
   perf <- performance(pred, "sens", "spec")
   sens <- perf@y.values
   spec <- perf@x.values
   cutoffs <- pred@cutoffs
   
   ppv <- npv <- NULL

   if(!is.na(prev)) {
      ppv <- sapply(1:length(cutoffs), function(j) {
         (sens[[j]] * prev) / (
            sens[[j]] * prev + (1 - spec[[j]]) * (1 - prev))
      })
      npv <- sapply(1:length(cutoffs), function(j) {
         (spec[[j]] * (1 - prev)) / (
            (1 - sens[[j]]) * prev + spec[[j]] * (1 - prev))
      })
   }

   list(
      pred=pred,
      perf=perf,
      sens=sens,
      spec=spec,
      cutoffs=cutoffs,
      ppv=ppv,
      npv=npv
   )
})

f <- sprintf("%s/%s.RData", outdir, name) 
cat("saving results to file", f, "\n")
save(res, file=f)


