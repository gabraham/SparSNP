#!/bin/bash

if [ -z "$1" ] ;
then
   echo "usage: auc.sh <numdirs> <title> [uni] [prev=%]"
   exit 1
fi

cat > .process.R <<EOF

library(ggplot2)

exper <- 1:$1
title <- "$2"
arg <- c("$3", "$4", "$5", "$6", "$7")
uni <- length(grep("^uni", arg)) >= 1
w <- grep("prev", arg)
prev <- if(length(w) > 0) {
   as.numeric(strsplit(arg[w], "=")[[1]][2])
} else NULL

w <- grep("minauc", arg)
minauc <- if(length(w) > 0) {
   as.numeric(strsplit(arg[w], "=")[[1]][2])
} else NULL

w <- grep("maxauc", arg)
maxauc <- if(length(w) > 0) {
   as.numeric(strsplit(arg[w], "=")[[1]][2])
} else NULL

w <- grep("minvar", arg)
minvar <- if(length(w) > 0) {
   as.numeric(strsplit(arg[w], "=")[[1]][2])
} else NULL

w <- grep("maxvar", arg)
maxvar <- if(length(w) > 0) {
   as.numeric(strsplit(arg[w], "=")[[1]][2])
} else NULL

w <- grep("minn", arg)
minn <- if(length(w) > 0) {
   as.numeric(strsplit(arg[w], "=")[[1]][2])
} else 0

w <- grep("maxn", arg)
maxn <- if(length(w) > 0) {
   as.numeric(strsplit(arg[w], "=")[[1]][2])
} else Inf

cat("prevalence:", prev, "\n")
cat("minauc:", minauc, "maxauc:", maxauc, "\n")
cat("minvar:", minvar, "maxvar:", maxvar, "\n")
cat("minn:", minn, "maxn:", maxn, "\n")

for(i in exper)
{
   setwd(sprintf("crossval%s", i))
   source("~/WTCCC/code/evalpred.R")
   setwd("../")
}

res <- lapply(exper, function(i) {
   cat("crossval", i, "\n")
   load(sprintf("crossval%s/crossval.RData", i))
   cv.auc.d
})
cv <- do.call(rbind, res)
cv <- cv[cv\$NonZero >= minn & cv\$NonZero <= maxn, ]

res.uni <- cv.uni <- NULL
if(uni)
{
   res.uni <- lapply(exper, function(i) {
      load(sprintf("crossval%s/crossval.RData", i))
      cv.uni.auc.d
   })
   cv.uni <- do.call(rbind, res.uni)
   cv.uni <- cv.uni[cv.uni\$NonZero >= minn & cv.uni\$NonZero <= maxn, ]
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
   source("~/WTCCC/code/varexp.R")
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

mytheme <- function(base_size=10)
{
   structure(list(
	 axis.text.x=theme_text(size=15),
	 axis.text.y=theme_text(size=20, hjust=1),
	 axis.title.x=theme_text(size=20),
	 axis.title.y=theme_text(size=20, angle=90),
	 plot.title=theme_text(size=30),
	 #axis.ticks=theme_blank(),
	 plot.title=theme_text(size=30),
	 legend.text=theme_text(size=20),
	 legend.title=theme_text(size=20, hjust=0),
	 legend.key.size=unit(2, "lines"),
	 legend.background=theme_rect(col=0, fill=0),
	 legend.key=theme_blank()
   ), class="options")
}

g <- if(uni) {
   ggplot(d, aes(x=NonZero, y=AUC, shape=Method, colour=Method))
} else {
   ggplot(d, aes(x=NonZero, y=AUC))
}

m <- round(max(log2(d\$NonZero)))
br <- 2^(0:m)
g <- g + geom_point(size=2.5)
g <- g + scale_x_log2("Number of SNPs in model", breaks=br, labels=br) 
g <- g + theme_bw() + mytheme()
g <- g + scale_colour_grey(start=0, end=0.5)
g <- g + stat_smooth(method="loess")
g <- g + coord_cartesian(ylim=c(minauc, maxauc))

pdf(sprintf("%s_AUC.pdf", title), width=14)
print(g)
dev.off()


if(!is.null(prev))
{
   g <- if(uni) {
      ggplot(d, aes(x=NonZero, y=VarExp, shape=Method, colour=Method))
   } else {
      ggplot(d, aes(x=NonZero, y=VarExp))
   }
   g <- g + geom_point(size=2.5) 
   g <- g + scale_x_log2("Number of SNPs in model", breaks=br, labels=br)
   g <- g + theme_bw() + mytheme()
   g <- g + scale_colour_grey(start=0, end=0.5)
   g <- g + stat_smooth(method="loess")
   g <- g + coord_cartesian(ylim=c(minvar, maxvar))
   
   
   pdf(sprintf("%s_VarExp.pdf", title), width=14)
   print(g)
   dev.off()

   m <- melt(d, measure.vars=c("AUC", "VarExp"))
   colnames(m)[colnames(m) == "variable"] <- "Measure"

   g <- if(uni) {
      ggplot(m, aes(x=NonZero, y=value, shape=Method, colour=Measure))
   } else {
      ggplot(m, aes(x=NonZero, y=value, colour=Measure))
   }
   g <- g + geom_point(size=2.5) 
   g <- g + scale_x_log2("Number of SNPs with non-zero weights",
      breaks=br, labels=br)
   
   g <- g + theme_bw() + mytheme()
   g <- g + scale_colour_grey(start=0, end=0.5)
   g <- g + stat_smooth(method="loess", size=2)
   g <- g + coord_cartesian(ylim=c(minvar, maxvar))
   g <- g + scale_y_continuous("Measure")
   g <- g + opts(legend.position="none")
   g <- g + geom_text(aes(x=512, y=0.5, label="VarExp"), size=10,
      colour="darkslategray")
   g <- g + geom_text(aes(x=512, y=0.75, label="AUC"), size=10, colour=1)
   g <- g + opts(plot.margin=unit(c(0, 0, 0.1, 1.5), "lines"))
   
   pdf(sprintf("%s_AUC_VarExp.pdf", title), width=8, height=6)
   print(g)
   grid.text("B", x=unit(0.025, "npc"), y=unit(0.94, "npc"),
	 gp=gpar(fontsize=50))
   dev.off()

   
}
EOF

Rscript .process.R
#/bin/rm .process.R

