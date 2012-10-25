#!/usr/bin/env Rscript

#options(error=NULL)

usage <- paste("usage: eval.R [title=<title>]",
      "[prev=<K>] [h2l=<V>] [mode=discovery|validation]",
      "(prev must be specified if h2l is specified)")

if(!exists("uni")){
   uni <- FALSE
}
best.uni.k <- best.uni <- best <- NULL

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
library(abind)

scale_x_log2 <- function(...)
{
   scale_x_continuous(..., trans=scales::log2_trans())
}

# Process prediction results

auc <- function(y, p)
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
   
   if(!exists("uni"))
      uni <- FALSE
   
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
	 y <- sapply(y, function(r) {
	    as.numeric(as.character(factor(y, labels=c(0, 1))))
	 })
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

      #browser()

      list(
         r=cbind(
	    lambda=lambda[1:nrow(res)],
	    Measure=res,
	    NonZero=nz[1:nrow(res), ]
         )#,
         #prediction=pred
      )
   })
   
   #pred <- lapply(cv, function(x) x$prediction)
   cv.d <- data.frame(do.call(rbind, lapply(cv, function(x) x$r)))
   
   # remove zeros
   #cv.d <- cv.d[cv.d$NonZero > 0, ] 

   cv.uni.d <- pred.uni <- NULL
   #if(uni)
   #{
   #   cv.uni <- lapply(folds, function(fold) {
   #      # should be 0/1
   #      y <- scan(sprintf("multivar_y.%02d", fold), quiet=TRUE)
   #   
   #      files <- sort(
   #         list.files(
   #	    pattern=sprintf(
   #	       "^multivar_beta.csv.[[:digit:]]+.%02d.pred[\\.bz2]*$", fold))
   #      )
   #   
   #      pr <- sapply(files, scan, quiet=TRUE)
   #   
   #      if(type == "AUC") {
   #         y <- as.numeric(as.character(factor(y, labels=c(0, 1))))
   #         res <- apply(pr, 2, auc, y=y)
   #      } else {
   #         res <- apply(pr, 2, r2, y=y)
   #      }

   #      #pred <- list(
   #      #   p=lapply(1:ncol(pr), function(j) pr[,j]),
   #      #   y=lapply(1:ncol(pr), function(j) y)
   #      #)
   #   
   #      # Multivar nonzero doesn't have intercept
   #      nz <- scan(sprintf("multivar_nonzero.csv.%02d", fold), quiet=TRUE)
   #      nz.thin <- scan(sprintf("multivar_nonzero.csv_thinned.%02d", fold),
   #         quiet=TRUE)
   #      status <- scan(
   #            list.files(pattern=sprintf("multivar_[irls|status]+.%02d",
   #     	     fold)), quiet=TRUE
   #      )

   #      w <- which(status == 3)
   #      if(length(w) == 0) {
   #         w <- length(status)
   #      }

   #      ok <- which(status[1:w] == 1)
   #
   #      list(
   #         r=cbind(
   #            Measure=res[ok],
   #            NonZeroPreThin=nz[ok],
   #            NonZero=nz.thin[ok]
   #         )#,
   #         #prediction=pred
   #      )
   #   })
   #   
   #   pred.uni <- lapply(cv.uni, function(x) x$prediction)
   #   cv.uni.d <- data.frame(
   #      do.call(rbind, lapply(cv.uni, function(x) x$r))
   #   )
   #}

   save(cv.d, cv.uni.d, 
      #pred, pred.uni,
      file="crossval.RData")
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
   pr <- abind(pr, along=3)
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

   d <- data.frame(NonZero=nz, Measure=res)
   d <- d[d$NonZero > 0, ]
   cv.d <- d

   ############################################################################    
   cv.uni.d <- NULL
   ## univariable/multivariable
   #files <- sort(
   #   list.files(
   #      pattern="^multivar_beta\\.csv\\.[[:digit:]]+\\.[[:digit:]]+\\.pred[\\.bz2]*$")
   #)
   # 
   ## remove predictions from bed models that didn't converge etc
   #lf <- list.files(pattern="multivar_[irls|status]+\\.[[:digit:]]+$",
   #      path=discovery.dir, full.names=TRUE)
   #status <- lapply(lf, scan, quiet=TRUE)
   #names(status) <- lf
   #f <- sapply(strsplit(lf, "\\."), tail, n=1)
   #l <- lapply(f, function(g) {
   #   grep(sprintf("^multivar_beta\\.csv\\.[[:digit:]]+\\.%s\\.pred$", g),
   #      files, value=TRUE)
   #})

   #l2 <- lapply(seq(along=l), function(i) {
   #   l[[i]][status[[i]] == 1]
   #})

   #l3 <- unlist(l2)

   ### not a sparse model
   #pr <- sapply(l3, scan, quiet=TRUE)
   #beta <- gsub("\\.pred", "", l3)
   ##nz <- sapply(beta, function(x) {
   ##   b <- scan(sprintf("%s/%s", discovery.dir, x), quiet=TRUE)
   ##   sum(b[-1] != 0)
   ##})
   #nz <- sapply(beta, function(x) {
   #   nrow(read.table(sprintf("%s/%s", discovery.dir, x))) - 1
   #})

   #if(type == "AUC") {
   #   y <- as.numeric(as.character(factor(y, labels=c(0, 1))))
   #   res <- apply(pr, 2, auc, y=y)
   #} else {
   #   res <- apply(pr, 2, r2, y=y)
   #}

   #d <- data.frame(NonZero=nz, Measure=res)
   #cv.uni.d <- d[d$NonZero > 0, ]

   save(cv.d, cv.uni.d, file="crossval.RData")
   list(cv.d, cv.uni.d)
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

params <- scan(sprintf("%s/discovery/params.txt", rootdir), what="character",
   quiet=TRUE)

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
#cv <- cv[cv$NonZero > 0, ]

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

#d <- cv[, vars]
d <- cv
#if(uni)
#{
#   cv.uni$Method <- "logistic"
#   if(measure == "AUC")
#   {
#      if(!is.null(prev)) {
#	 cv.uni$VarExp <- varexp(K=prev,
#	    auc=cv.uni$Measure)[, "varexp"]
#      }
#
#      if(!is.null(h2l)) {
#	 cv.uni$GenVarExp <- varexp(K=prev,
#	    auc=cv.uni$Measure, h2l=h2l)[, "genvarexp"]
#      }
#   }
#
#   #d <- rbind(d, cv.uni[cv.uni$Status == 1, vars])
#   d <- rbind(d, cv.uni[, vars])
#}

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

numtasks <- length(grep("^NonZero_[[:digit:]]+", colnames(d)))

#if(numtasks > 0) {
#   d2 <- lapply(1:numtasks, function(k) {
#      v1 <- grep(sprintf("%s_%s$", k), measure, colnames(d), value=TRUE)
#      v2 <- grep(sprintf("NonZero_%s$", k), colnames(d), value=TRUE)
#      data.frame(d$lambda, Measure=d[, v1], NonZero=d[,v2], d$Method, Task=k)
#   })
#   d3 <- do.call(rbind, d2)
#} else {
#   d3 <- data.frame(d, Task=1)
#}
#
#d3 <- d3[d3$NonZero > 0 ,]
#d3$Task <- factor(d3$Task)
#
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

g <- if(uni) {
   ggplot(d, aes(x=NonZero, y=Measure, shape=Method, colour=Method))
} else {
   ggplot(d, aes(x=NonZero, y=Measure))
}

m <- round(max(log2(d[, grep("^NonZero", colnames(d))])))
br <- 2^(0:m)
br <- br[br <= max(d$NonZero)]
expr <- ifelse(measure == "AUC", "AUC", expression(R^2))
g <- g + geom_point(size=2.5)
#g <- g + scale_x_log2("Number of SNPs in model", breaks=br, labels=br) 
g <- g + scale_x_log2("Number of SNPs in model") 
g <- g + scale_y_continuous(expr)
g <- g + theme_bw() + mytheme
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
   g <- g + theme_bw() + mytheme
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
   g <- g + theme_bw() + mytheme
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
         s <- scan(sprintf("crossval%s/nonzero.csv.%02d", rep, fold - 1),
	    quiet=TRUE)
	 
	 # Never select zero as a legitimate model
	 s[s == 0] <- -Inf
         w <- which.min((s - m)^2)
         #cat("best:", s[w], " ")
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

topsnps <- NULL

if(mode == "discovery")
{
   res <- tabulate.snps(best, cv)
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
   write.table(topsnps, file="topsnps.txt", quote=FALSE, row.names=FALSE)

   if(uni)
   {
      res.uni <- tabulate.snps(best.uni, cv.uni)
      snps.uni <- res.uni$snps
      best.uni <- res.uni$best
      best.uni.k <- res.uni$k
      
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

save(cv, topsnps, best, best.uni, best.k, best.uni.k,
   file=sprintf("%s.RData", title))

