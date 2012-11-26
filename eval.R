#!/usr/bin/env Rscript

#options(error=NULL)

usage <- paste("usage: eval.R [title=<title>]",
      "[prev=<K>] [h2l=<V>] [mode=discovery|validation]",
      "(prev must be specified if h2l is specified)")

args <- commandArgs(TRUE)
s <- strsplit(args, "=")
for(m in s) {
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

#if(!is.null(best))
#   best <- as.integer(best)
best <- NULL

library(ggplot2)
library(scales)
library(grid)
library(abind)

scale_x_log2 <- function(...)
{
   scale_x_continuous(..., trans=scales::log2_trans())
}

# Process prediction results

auc <- function(p, y)
{
   y <- cbind(y)
   p <- cbind(p)

   sapply(1:ncol(y), function(k) {
      r <- rank(p[,k])
      n1 <- sum(y[,k] == 1)
      n0 <- length(y[,k]) - n1
      u <- sum(r[y[,k] == 1]) - n1 * (n1 + 1) / 2
      u / (n0 * n1)
   })
}

# R^2
#
# Same definition as in summary.lm:
#  R^2 = 1 - Sum(R[i]^2) / Sum((y[i]- y*)^2),
#  where R are the residuals and y* is the intercept 
#
# f: the fitted (predicted) value of y
# y: the actual value of y
r2 <- function(f, y)
{
   f <- cbind(f)
   y <- cbind(y)

   sapply(1:ncol(y), function(k) {
      1 - sum((f[,k] - y[,k])^2) / sum((y[,k] - mean(y[,k]))^2)
   })
}

# for cross-validation
evalpred.crossval <- function(type=NULL, dir=NULL)
{
   type <- match.arg(type, c("AUC", "R2"))

   folds <- scan("folds.txt", quiet=TRUE)
   folds <- unique(sort(folds))
   
   cv.d <- pred <- NULL
   cv <- lapply(folds, function(fold) {
      y <- as.matrix(read.table(sprintf("y.%02d", fold), sep=","))
      
      if(length(unique(y)) == 1) {
	 stop("only one unique value of y found, cannot evaluate prediction",
	    "(file=", sprintf("y.%02d", fold), ")")
      }
   
      files <- sort(
         list.files(
	    pattern=sprintf(
	       "^beta\\.csv\\.[[:digit:]]+\\.%02d\\.pred[\\.bz2]*$", fold))
      )
   
      pr <- lapply(files, function(f) as.matrix(read.table(f, sep=",")))
      lambda <- scan(sprintf("lambda1path.csv.%02d", fold), quiet=TRUE)
      nz <- as.matrix(read.table(sprintf("nonzero.csv.%02d", fold), sep=","))

      pr <- abind(pr, along=3)
      
      # Don't count the intercept, except for first #nonzero
      #nz[-1] <- nz[-1] - 1  
   
      res <- if(type == "AUC") {
	 # y <- sapply(length(, function(r) {
	 #    as.numeric(as.character(factor(y, labels=c(0, 1))))
	 # })
	 y <- cbind(as.numeric(as.character(factor(y, labels=c(0, 1)))))
	 apply(pr, 3, function(p) {
	    auc(p, y)
	 })
      } else {
	 apply(pr, 3, function(p) {
	    r2(p, y)
	 })
      }
      
      # apply messes up orientation when res is a vector
      res <- t(rbind(res)) 

      #pred <- list(
      #   p=lapply(1:ncol(pr), function(j) pr[,j]),
      #   y=lapply(1:ncol(pr), function(j) y)
      #)

      if(ncol(y) > 1) {
	 colnames(res) <- paste("Measure", 1:ncol(y), sep="_")
	 colnames(nz) <- paste("NonZero", 1:ncol(y), sep="_")
      } else {
	 colnames(res) <- "Measure"
	 colnames(nz) <- "NonZero"
      }

      list(
         r=cbind(
	    lambda=lambda[1:nrow(res)],
	    Measure=res,
	    NonZero=nz[1:nrow(res), ]
         )#,
         #prediction=pred
      )
   })
   
   cv.d <- data.frame(do.call(rbind, lapply(cv, function(x) x$r)))
   
   save(cv.d, file="crossval.RData")
}

# for testing on validation dataset
evalpred.validation <- function(type=NULL, discovery.dir=NULL)
{
   type <- match.arg(type, c("AUC", "R2"))
   if(is.null(discovery.dir))
      stop("discovery.dir not specified")

   y <- scan("y.txt", quiet=TRUE)

   ###########################################################################    
   # lasso
   files <- sort(
      list.files(
         pattern="^beta\\.csv\\.[[:digit:]]+\\.[[:digit:]]+\\.pred[\\.bz2]*$")
   )
   
   pr <- lapply(files, function(f) as.matrix(read.table(f, sep=",")))
   pr <- abind(pr, along=2)
   beta <- gsub("\\.pred", "", files)
   nz <- sapply(beta, function(x) {
      nrow(read.table(sprintf("%s/%s", discovery.dir, x))) - 1
   })
   
   if(type == "AUC") {
      y <- as.numeric(as.character(factor(y, labels=c(0, 1))))
      res <- apply(pr, 2, auc, y=y)
   } else {
      res <- apply(pr, 2, r2, y=y)
   }

   cv.d <- data.frame(NonZero=nz, Measure=res)
   cv.d <- cv.d[cv.d$NonZero > 0, ]
 
   save(cv.d, file="crossval.RData")
   list(cv.d)
}


# Genetic variance explained
# Wray et al "Genetic interpretation of Area under the ROC Curve in Genomic
# Profiling", PLoS Genetics, 2:e1000864, 2010
#
#
# K: prevelance > 0
# h2l: heritability
# auc: Area under ROC Curve
# 
# Must specify K and at least one of h2l, auc
#
# Output:
#   K:      prevalence as supplied
#   h2l:    genetic heritability as supplied
#   aucmax: maximum AUC given K and h2l (explains all genetic variance)
#   auc:    AUC as supplied
#
varexp <- function(K=NULL, h2l=NULL, auc=NULL, auc.se=NULL)
{
   if(is.null(K) || all(is.null(c(h2l, auc)))) {
      stop("Must specify K and at least one of h2l, auc")
   }

   if(any(c(K, h2l, auc) <= 0) || any(c(K, h2l, auc) > 1)) {
      stop("K, h2l, and auc must be in the range (0, 1]")
   }

   T <- qnorm(1 - K)
   z <- dnorm(T)
   i <- z / K
   v <- -i * K / (1 - K)
   
   aucmax <- as.numeric(NA)
   if(!is.null(h2l)) {
      # maximum achievable AUC
      mean1 <- i * h2l
      mean2 <- v * h2l
      var1 <- h2l * (1 - h2l * i * (i - T))
      var2 <- h2l * (1 - h2l * v * (v - T))
      
      d <- (mean1 - mean2) / sqrt(var1 + var2)
      aucmax <- pnorm(d)
   }
   
   h2lx <- rho2gg <- h2lx.se <- as.numeric(NA)
   if(!is.null(auc)) {
      # variance of genetic profile
      Q <- qnorm(auc, lower.tail=TRUE)
      h2lx <- 2 * Q^2 / ((v - i)^2 + Q^2 * i * (i - T) + v * (v - T))

      if(!is.null(h2l)) {
	 rho2gg <- h2lx / h2l
	 rho2gg <- pmin(rho2gg, 1)
      }

      # Based on delta rule SE(varexp(AUC)) = SE(AUC) * varexp(AUC)'
      if(!is.null(auc.se)) {
	 dQ <- 1 / dnorm(Q)
      	 h2lxd <- -8 * Q^2 * dQ^2 * i * (i - T) / (
      	    ((v - i)^2 + Q^2 * i * (i - T) + v * (v - T))^2
      	 )
	 h2lx.se <- auc.se * h2lxd
      }
   }

   cbind(K=K, h2l=h2l, aucmax=aucmax, varexp=h2lx, genvarexp=rho2gg,
      varexp.se=h2lx.se)
}

cat("mode:", mode, "\n")

rootdir <- getwd()

if(mode == "discovery") {
   fun <- evalpred.crossval
   setwd("discovery")
} else {
   fun <- evalpred.validation
   setwd("validation")
}

params <- scan(sprintf("%s/discovery/params.txt", rootdir),
   what="character", quiet=TRUE)

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

d <- cv

mytheme <- theme(
	 axis.text.x=element_text(size=15),
	 axis.text.y=element_text(size=20, hjust=1),
	 axis.title.x=element_text(size=20),
	 axis.title.y=element_text(size=20, angle=90),
	 plot.title=element_text(size=30),
	 #axis.ticks=element_blank(),
	 plot.title=element_text(size=30),
	 legend.text=element_text(size=20),
	 legend.title=element_text(size=20, hjust=0),
	 legend.key.size=unit(2, "lines"),
	 legend.background=element_rect(colour=0, fill=0),
	 legend.key=element_blank()
)

numtasks <- max(length(grep("^NonZero_[[:digit:]]+", colnames(d))), 1)

if(numtasks > 1) {
   d2 <- lapply(1:numtasks, function(k) {
      v1 <- grep(sprintf("Measure_%s$", k), colnames(d), value=TRUE)
      v2 <- grep(sprintf("NonZero_%s$", k), colnames(d), value=TRUE)
      data.frame(lambda=d$lambda, Measure=d[, v1], NonZero=d[,v2], Method=d$Method, Task=k)
   })
   d3 <- do.call(rbind, d2)
} else {
   d3 <- data.frame(d, Task=1)
}

d3 <- d3[d3$NonZero > 0 ,]
d3$Task <- factor(d3$Task)

#g <- ggplot(d3, aes(x=NonZero, y=Measure, colour=Task, line=Task))
#g <- g + geom_smooth(method="loess", size=1.5) + scale_x_log10()   
#g <- g + theme_bw() + mytheme + scale_y_continuous(expression(R^2))
#g <- g + guides(colour=guide_legend(ncol=2))
##g <- g + geom_point(alpha=0.1)
#
#pdf(sprintf("%s.pdf", measure), width=14)
#print(g)
#dev.off()

## "Loss" pooled over all tasks
#if(numtasks > 0) {
#   g <- ggplot(d3, aes(x=NonZero, y=Measure))
#   g <- g + geom_smooth(method="loess", se=FALSE, size=1.5) + scale_x_log10()   
#   g <- g + theme_bw() + mytheme + scale_y_continuous(expression(R^2))
#   #g <- g + guides(colour=guide_legend(ncol=2))
#   
#   lo <- loess(Measure ~ log(NonZero), data=d3)
#   s <- predict(lo, newdata=data.frame(NonZero=1:max(d3$NonZero)))
#   m <- which(s == max(s, na.rm=TRUE))
#   cat("Best model at", m, "SNPs with predictive measure", s[m], "\n") 
#   
#   g <- g + geom_vline(xintercept=m, lty=2)
#   g <- g + annotate("text", x=m + 15, y=-0.05, label=m)
#   
#   pdf(sprintf("%s_all.pdf", measure), width=14)
#   print(g)
#   dev.off()
#}

g <- ggplot(d3, aes(x=NonZero, y=Measure, colour=Task))

m <- round(max(log2(d3$NonZero)))
br <- 2^(0:m)
br <- br[br <= max(d3$NonZero)]
expr <- ifelse(measure == "AUC", "AUC", expression(R^2))
#g <- g + geom_point(size=2.5, alpha=0.2)
#g <- g + scale_x_log2("Number of SNPs in model", breaks=br, labels=br) 
g <- g + scale_x_log2("Number of SNPs in model") 
g <- g + scale_y_continuous(expr)
g <- g + theme_bw() + mytheme
#g <- g + scale_colour_grey(start=0, end=0.5)
g <- g + stat_smooth(method="loess", size=2)

pdf(sprintf("%s_%s.pdf", title, measure), width=14)
print(g)
dev.off()


if(!is.null(prev))
{
   g <- ggplot(d, aes(x=NonZero, y=VarExp))
   g <- g + geom_point(size=2.5) 
   #g <- g + scale_x_log2("Number of SNPs in model", breaks=br, labels=br)
   g <- g + scale_x_log2("Number of SNPs in model")
   g <- g + theme_bw() + mytheme
   g <- g + scale_colour_grey(start=0, end=0.5)
   g <- g + stat_smooth(method="loess")
   
   pdf(sprintf("%s_VarExp.pdf", title), width=14)
   print(g)
   dev.off()
}

if(!is.null(h2l))
{
   g <- ggplot(d, aes(x=NonZero, y=GenVarExp))
   g <- g + geom_point(size=2.5) 
   #g <- g + scale_x_log2("Number of SNPs in model", breaks=br, labels=br)
   g <- g + scale_x_log2("Number of SNPs in model")
   g <- g + theme_bw() + mytheme
   g <- g + scale_colour_grey(start=0, end=0.5)
   g <- g + stat_smooth(method="loess")
   
   pdf(sprintf("%s_GenVarExp.pdf", title), width=14)
   print(g)
   dev.off()
}

save(d3, file="discovery.RData")

if(numtasks > 1)
   stop("SNP tabulation broken for multitask")

# Find SNPs for best model in cross-validation
# best: if specified as an integer >0 , then return the SNPs at size model
#       size. Otherwise, best model will be automatically determined using the
#       highest "Measure" (AUC or R^2)
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

	 s <- as.matrix(read.table(
	    sprintf("crossval%s/nonzero.csv.%02d", rep, fold - 1),
	    sep=","))
	 
	 # Never select zero as a legitimate model
	 s[s == 0] <- -Inf
         w <- which.min((s - m)^2)
         b <- read.table(
	    sprintf("crossval%s/beta.csv.%02d.%02d",
	       rep, w - 1, fold - 1), sep=",")
	 # exclude intercept
         b[-1, 1]
      })
   })

   b <- ifelse(is.null(best), s[m], best)

   list(best=b, snps=sort(table(unlist(l)), decreasing=TRUE), k=s[m])
}

get_topsnps <- function(...)
{
   res <- tabulate.snps(...)
   snps <- res$snps
   best <- res$best
   best.k <- res$k
   
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
   topsnps
}

if(mode == "discovery")
{
   topsnps <- if(numtasks == 1) {
      s <- get_topsnps(best, d3)
      write.table(s, file="topsnps.txt", quote=FALSE, row.names=FALSE)
      s
   } else {
     lapply(1:numtasks, function(k) {
	 s <- get_topsnps(best, d3[d3$Task == k, ])
         write.table(s, file=sprintf("topsnps_%d.txt", k),
	    quote=FALSE, row.names=FALSE)
	 s
      })
   }
}

# Change from generic name to actual name (AUC/R2)
colnames(cv)[colnames(cv) == "Measure"] <- measure

cv <- d

save(cv, d3, topsnps, #best, best.k, 
   file=sprintf("%s.RData", title))

