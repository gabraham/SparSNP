#!/bin/bash

if [ -z "$1" ] ;
then
   echo "usage: auc.sh <numdirs> <title> [uni] [prev=%]"
   exit 1
fi

cat > .process.R <<EOF

library(ggplot2)
library(Hmisc)

exper <- 1:$1
title <- "$2"
arg <- c("$3", "$4")
uni <- length(grep("^uni", arg)) >= 1
w <- grep("prev", arg)
prev <- if(length(w) > 0) {
   as.numeric(strsplit(arg[w], "=")[[1]][2])
} else NULL

cat("prevalence:", prev, "\n")


for(i in exper)
{
   setwd(sprintf("crossval%s", i))
   source("../evalpred.R")
   setwd("../")
}

res <- lapply(exper, function(i) {
   cat("crossval", i, "\n")
   load(sprintf("crossval%s/crossval.RData", i))
   cv.auc.d
})
cv <- do.call(rbind, res)

res.uni <- cv.uni <- NULL
if(uni)
{
   res.uni <- lapply(exper, function(i) {
      load(sprintf("crossval%s/crossval.RData", i))
      cv.uni.auc.d
   })
   cv.uni <- do.call(rbind, res.uni)
}

if(title == "")
   title <= "combined.pdf"

save(res, cv, res.uni, cv.uni, file=sprintf("%s.RData", title))

cv\$Method <- "lasso"

vars <- if(is.null(prev)) {
   c("NonZero", "AUC", "Method")
} else {
   c("NonZero", "AUC", "Method", "VarExp")
}

if(!is.null(prev))
{
   source("varexp.R")
   cv\$VarExp <- varexp(K=prev, auc=cv\$AUC)[, "varexp"]
   if(uni) {
      cv.uni\$VarExp <- varexp(K=prev, auc=cv.uni\$AUC)[, "varexp"]
   }
}

d <- cv[, vars]
if(uni)
{
   cv.uni\$Method <- "logistic"
   d <- rbind(d, cv.uni[cv.uni\$Status == 1, vars])
}


g <- ggplot(d, aes(x=NonZero, y=AUC, shape=Method, colour=Method))
g <- g + geom_point(size=2.5) + scale_x_log2() + theme_bw()
g <- g + scale_colour_grey(start=0, end=0.5)
g <- g + stat_smooth(method="loess")

pdf(sprintf("%s_AUC.pdf", title), width=12)
print(g)
dev.off()


if(!is.null(prev))
{
   g <- ggplot(d, aes(x=NonZero, y=VarExp, shape=Method, colour=Method))
   g <- g + geom_point(size=2.5) + scale_x_log2() + theme_bw()
   g <- g + scale_colour_grey(start=0, end=0.5)
   g <- g + stat_smooth(method="loess")
   
   pdf(sprintf("%s_VarExp.pdf", title), width=12)
   print(g)
   dev.off()
}
EOF

Rscript .process.R
#/bin/rm .process.R

