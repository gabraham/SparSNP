#!/usr/bin/env Rscript

options(error=NULL)

usage <- paste("usage: eval.R [title=<title>]",
      "[prev=<K>] [h2l=<V>] [mode=discovery|validation]",
      "(prev must be specified if h2l is specified)")

if(!exists("uni")){
   uni <- FALSE
}
best.uni <- best <- NULL

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
   if(is.na(prev) || prev < 0 || prev > 1.0) {
      stop(usage)
   }
} else {
   prev <- NULL
}

if(exists("h2l")) {
   if(is.null(prev)) {
      stop(usage)
   }
   h2l <- as.numeric(h2l)
   if(is.na(h2l) || h2l < 0 || h2l > 1.0)
      stop(usage)
} else {
   h2l <- NULL
}

if(!is.null(best))
   best <- as.integer(best)
if(!is.null(best.uni))
   best.uni <- as.integer(best.uni)

library(ggplot2)
library(scales)
library(grid)

scale_x_log2 <- function(...)
{
   scale_x_continuous(..., trans=scales::log2_trans())
}

source("~/Code/cd/src/sparsnp/evalpred.R")
source("~/Code/cd/src/sparsnp/varexp.R")


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

extract.params <- function(nm)
{  
   strsplit(params[grep(nm, params)], "=")[[1]][2]
}

model <- extract.params("MODEL")
nfolds <- as.integer(extract.params("NFOLDS"))
nreps <- as.integer(extract.params("NREPS"))

measure <- if(model == "sqrhinge" || model == "logistic") {
   "AUC"
} else if(model == "linear") {
   "R2"
} else stop("unknown model in params.txt:", model)

if(measure != "AUC" && !is.null(prev))
{
   prev <- NULL
   h2l <- NULL
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


cv$Method <- "lasso"

vars <- if(is.null(prev)) {
   c("NonZero", "Measure", "Method")
} else if(is.null(h2l)) {
   c("NonZero", "Measure", "Method", "VarExp")
} else {
   c("NonZero", "Measure", "Method", "VarExp", "GenVarExp")
}

if(measure == "AUC")
{
   if(!is.null(prev)) {
      cv$VarExp <- varexp(K=prev, auc=cv$Measure)[, "varexp"]
   }

   if(!is.null(h2l)) {
      cv$GenVarExp <- varexp(K=prev, auc=cv$Measure, h2l=h2l)[, "genvarexp"]
   }
}

d <- cv[, vars]
if(uni)
{
   cv.uni$Method <- "logistic"
   if(measure == "AUC")
   {
      if(!is.null(prev)) {
	 cv.uni$VarExp <- varexp(K=prev,
	    auc=cv.uni$Measure)[, "varexp"]
      }

      if(!is.null(h2l)) {
	 cv.uni$GenVarExp <- varexp(K=prev,
	    auc=cv.uni$Measure, h2l=h2l)[, "genvarexp"]
      }
   }

   #d <- rbind(d, cv.uni[cv.uni$Status == 1, vars])
   d <- rbind(d, cv.uni[, vars])
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
br <- br[br <= max(d$NonZero)]
g <- g + geom_point(size=2.5)
#g <- g + scale_x_log2("Number of SNPs in model", breaks=br, labels=br) 
g <- g + scale_x_log2("Number of SNPs in model") 
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
   #g <- g + scale_x_log2("Number of SNPs in model", breaks=br, labels=br)
   g <- g + scale_x_log2("Number of SNPs in model")
   g <- g + theme_bw() + mytheme()
   g <- g + scale_colour_grey(start=0, end=0.5)
   g <- g + stat_smooth(method="loess")
   
   pdf(sprintf("%s_VarExp.pdf", title), width=14)
   print(g)
   dev.off()
}

if(!is.null(h2l))
{
   g <- if(uni) {
      ggplot(d, aes(x=NonZero, y=GenVarExp, shape=Method, colour=Method))
   } else {
      ggplot(d, aes(x=NonZero, y=GenVarExp))
   }
   g <- g + geom_point(size=2.5) 
   #g <- g + scale_x_log2("Number of SNPs in model", breaks=br, labels=br)
   g <- g + scale_x_log2("Number of SNPs in model")
   g <- g + theme_bw() + mytheme()
   g <- g + scale_colour_grey(start=0, end=0.5)
   g <- g + stat_smooth(method="loess")
   
   pdf(sprintf("%s_GenVarExp.pdf", title), width=14)
   print(g)
   dev.off()
}


# Find SNPs for best model in cross-validation
tabulate.snps <- function(best=NULL, d)
{
   if(is.null(best))
   {
      lo <- loess(Measure ~ log(NonZero), data=d)
      s <- predict(lo, newdata=data.frame(NonZero=1:max(d$NonZero)))
      m <- which(s == max(s, na.rm=TRUE))
      cat("Best model at", m, "SNPs with predictive measure", s[m], "\n") 
   }
   else
   {
      m <- best
      cat("tabulate.snps: best=", best, "\n")
   }

   l <- lapply(1:nreps, function(rep) {
      l <- lapply(1:nfolds, function(fold) {
         s <- scan(sprintf("crossval%s/nonzero.csv.%02d", rep, fold - 1))
	 
	 # Never select zero as a legitimate model
	 s[s == 0] <- -Inf
         w <- which.min((s - m)^2)
         cat("best:", s[w], "\n")
         b <- read.table(
	    sprintf("crossval%s/beta.csv.%02d.%02d",
	       rep, w - 1, fold - 1), sep=":")
	 # exclude intercept
         b[-1, 1]
      })
   })

   b <- ifelse(is.null(best), s[m], best)
   
   list(best=b, snps=sort(table(unlist(l)), decreasing=TRUE))
}

topsnps <- NULL

if(mode == "discovery")
{
   res <- tabulate.snps(best, cv)
   snps <- res$snps
   best <- res$best
   
   topsnps <- cbind("NA"=numeric(0))
   
   if(length(snps) > 0)
   {
      #rs <- scan("snps.txt", what=character())
      ref <- read.table("snps.txt", header=FALSE, sep="",
	    stringsAsFactors=FALSE)
      rs <- ref[,1]
      names(snps) <- rs[as.integer(names(snps))]
      
      topsnps <- data.frame(
         RS=rownames(snps),
         Counts=snps,
         Proportion=snps / nreps / nfolds,
         Replications=nreps * nfolds
      )
   }
   write.table(topsnps, file="topsnps.txt", quote=FALSE, row.names=FALSE)

   if(uni)
   {
      res.uni <- tabulate.snps(best.uni, cv.uni)
      snps.uni <- res.uni$snps
      best.uni <- res.uni$best
      
      topsnps.uni <- cbind("NA"=numeric(0))
      
      if(length(snps.uni) > 0)
      {
         #rs.uni <- scan("snps.txt", what=character())
	 ref <- read.table("snps.txt", header=FALSE, sep="",
	       stringsAsFactors=FALSE)
	 rs.uni <- ref[,1]
         names(snps.uni) <- rs[as.integer(names(snps.uni))]
         
         topsnps.uni <- data.frame(
            RS=rownames(snps.uni),
            Counts=snps.uni,
            Proportion=snps.uni / nreps / nfolds,
            Replications=nreps * nfolds
         )
      }
      write.table(topsnps.uni, file="topsnps_uni.txt",
	 quote=FALSE, row.names=FALSE)
   }
}

# Change from generic name to actual name (AUC/R2)
colnames(cv)[colnames(cv) == "Measure"] <- measure

cv <- d

save(cv, topsnps, best, best.uni, file=sprintf("%s.RData", title))

