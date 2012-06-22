#!/usr/bin/env Rscript

usage <- paste("usage: evalprofile.R",
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

#target <- r[1]
#target.base <- r[2]
#famf <- r[3]
#prev <- as.numeric(r[4])
#dir <- r[5]

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
   
   score <- dither(prof$SCORE + intercept)
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

#save(pheno, score, intercept, prev, prof, pred, perf, sens, spec,
#   ppv, npv, file=sprintf("%s.RData", target.base))
f <- sprintf("%s/%s.RData", outdir, name) 
cat("saving results to file", f, "\n")
save(res, file=f)


