#!/usr/bin/env Rscript

usage <- "usage: eval.R [title=<title>] [prev=<K>] [mode=discovery|validation]"

uni <- FALSE

args <- commandArgs(TRUE)
s <- strsplit(args, "=")
for(m in s)
{
   eval(parse(text=sprintf("%s=\"%s\"", m[1], m[2])))
}

if(!exists("mode", mode="character"))
   mode <- "discovery"

if(!mode %in% c("discovery", "validation"))
   stop(usage)

if(!exists("title", mode="character"))
   title <- mode

if(exists("prev")) {
   prev <- as.numeric(prev)
   if(is.na(prev) || prev < 0 || prev > 1.0)
      stop(usage)
} else {
   prev <- NULL
}

library(ggplot2)
source("evalpred.R")
source("varexp.R")


cat("mode:", mode, "\n")

rootdir <- getwd()

if(mode == "discovery") {
   fun <- evalpred.crossval
   setwd("discovery")
} else {
   fun <- evalpred.validation
   setwd("validation")
}

params <- scan(sprintf("%s/discovery/params.txt", rootdir), what="character")
model <- strsplit(params[grep("MODEL", params)], "=")[[1]][2]

measure <- if(model == "sqrhinge") {
   "AUC"
} else if(model == "linear") {
   "R2"
} else stop("unknown model in params.txt:", model)

if(measure != "AUC" && !is.null(prev))
{
   prev <- NULL
   cat("Model is", model, "; ignoring prev parameter\n")
}

dirs <- list.files(pattern="^crossval[[:digit:]]+$")
for(d in dirs)
{
   setwd(d)
   fun(type=measure, sprintf("%s/discovery/%s", rootdir, d))
   setwd("../")
}

res <- lapply(seq(along=dirs), function(i) {
   cat("crossval", i, "\n")
   load(sprintf("crossval%s/crossval.RData", i))
   cv.d
})
cv <- do.call(rbind, res)

res.uni <- cv.uni <- NULL
if(uni)
{
   res.uni <- lapply(seq(along=dirs), function(i) {
      load(sprintf("crossval%s/crossval.RData", i))
      cv.uni.d
   })
   cv.uni <- do.call(rbind, res.uni)
}

save(res, cv, file=sprintf("%s.RData", title))

cv$Method <- "lasso"

vars <- if(is.null(prev)) {
   c("NonZero", "Measure", "Method")
} else {
   c("NonZero", "Measure", "Method", "VarExp")
}

if(!is.null(prev) && measure == "AUC")
{
   cv$VarExp <- varexp(K=prev, auc=cv$Measure)[, "varexp"]
}

d <- cv[, vars]
if(uni)
{
   cv.uni$Method <- "logistic"
   d <- rbind(d, cv.uni[cv.uni$Status == 1, vars])
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
   ggplot(d, aes(x=NonZero, y=Measure, shape=Method, colour=Method))
} else {
   ggplot(d, aes(x=NonZero, y=Measure))
}

m <- round(max(log2(d$NonZero)))
br <- 2^(0:m)
g <- g + geom_point(size=2.5)
g <- g + scale_x_log2("Number of SNPs in model", breaks=br, labels=br) 
g <- g + scale_y_continuous(measure)
g <- g + theme_bw() + mytheme()
g <- g + scale_colour_grey(start=0, end=0.5)
g <- g + stat_smooth(method="loess")

pdf(sprintf("%s_%s.pdf", title, measure), width=14)
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
   
   pdf(sprintf("%s_VarExp.pdf", title), width=14)
   print(g)
   dev.off()
}

