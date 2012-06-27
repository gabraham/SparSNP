#!/usr/bin/env Rscript

# Fits an unpenalised model to the SNPs selected by the lasso


args <- commandArgs(TRUE)
s <- strsplit(args, "=")
for(m in s)
{
   eval(parse(text=sprintf("%s=\"%s\"", m[1], m[2])))
}

if(!exists("probe"))
{
   stop("secondstage.R: must specify probe=XXXX")
}

probe <- as.character(probe)

load("~/Data/DILGOM/common.RData")

if(!probe %in% colnames(expr.regress))
{
   stop("probe ", probe, " not found in data")
}

library(MASS)

lf <- list.files(
   pattern="^beta\\.csv\\.[[:digit:]]+\\.[[:digit:]]+[\\.bz2]*$",
   path=".")

v <- lapply(lf, function(f) {
   read.table(f, sep=":")[-1, 1]
})
nz <- sapply(v, length)

vt <- as.integer(names(table(unlist(v))))

root <- "~/Data/DILGOM/Human610/genotypes_qc_imputed_common_human610"
bim <- read.table(sprintf("%s.bim", root), header=FALSE,
      stringsAsFactors=FALSE)
rownames(bim) <- bim[, 2]

write.table(bim[vt, 2], file="snps.txt", col.names=FALSE, row.names=FALSE,
      quote=FALSE)

system(paste(
   "plink --noweb", "--bfile", root, "--recodeA",
   "--extract snps.txt"
))

d <- read.table("plink.raw", header=TRUE, stringsAsFactors=FALSE)
nm.dilgom <- colnames(d)[-(1:6)]
x <- as.matrix(d[, -(1:6)])
w <- is.na(x)
x[w] <- sample(0:2, sum(w), TRUE)
colnames(x) <- gsub("_[AGCT]", "", colnames(x))
#y <- expr.regress[, probe]

#fam <- read.table(sprintf("%s.fam", root), header=FALSE,
#      stringsAsFactors=FALSE)
fam <- read.table(
   "~/Data/DILGOM/Human610/Test/genotypes_qc_imputed_common_human610_ILMN_1688423.fam",
   header=FALSE, stringsAsFactors=FALSE)
y <- fam[,6]
#y <- rnorm(nrow(fam))

system(paste(
   "plink --noweb --bfile ~/Data/DILGOM/Human610/Test/test",
   "--extract snps.txt --recodeA --out test"
))
d.val <- read.table("test.raw", header=TRUE, stringsAsFactors=FALSE)
x.val <- as.matrix(d.val[, -(1:6)])
w <- is.na(x.val)
x.val[w] <- sample(0:2, sum(w), TRUE)
colnames(x.val) <- gsub("_[AGCT]", "", colnames(x.val))
y.val <- d.val[, 6]
#b <- drop(ginv(cbind(1, x[, s]))) %*% y
#names(b) <- c("(Intercept)", s)
#p.test <- drop(cbind(1, x.test[, s]) %*% b)
#r2 <- 1 - sum((p.test - y.test)^2) / sum((y.test - mean(y.test))^2)

folds <- scan("folds.txt")

res <- lapply(1:max(folds), function(fold) {
   xtrain <- x[folds != fold, ]
   ytrain <- y[folds != fold]
   xtest <- x[folds == fold, ]
   ytest <- y[folds == fold]

   lf <- list.files(
      pattern=sprintf("^beta\\.csv\\.[[:digit:]]+\\.%02d+[\\.bz2]*$", fold - 1),
      path=".")
   res <- lapply(lf, function(f) {
      v <- read.table(f, sep=":")[-1, 1]
      nz <- length(v)
      if(length(v) == 0) {
	 c(NonZero=nz, R2=0, R2val=0)
      } else {
	 rs <- bim[v, 2]
	 b <- ginv(cbind(1, xtrain[, rs])) %*% ytrain

	 p <- cbind(1, xtest[, rs]) %*% b
	 r2 <- 1 - sum((p - ytest)^2) / sum((ytest - mean(ytest))^2)

	 p.val <- cbind(1, x.val[, rs]) %*% b
	 r2.val <- 1 - sum((p.val - y.val)^2) / sum((y.val - mean(y.val))^2)
	 
	 c(NonZero=nz, R2=r2, R2val=r2.val)
      }
   })
   do.call(rbind, res)
})

res <- data.frame(do.call(rbind, res))
stop()

w <- which(nz >= 150 & nz <= 250)
vs <- sort(table(unlist(v[w])), decreasing=TRUE)
rs <- vs / length(w)
names(rs) <- bim[as.integer(names(vs)), 2]

pdf("discovery_2nd_stage.pdf")
plot(R2 ~ NonZero, data=res, main=probe)
abline(h=0, lty=2)
dev.off()


#egc <- function()
#{
#   load("~/Data/EGC/ForGadEGC/common.RData")
#   scale(expr[, probe], center=TRUE, scale=FALSE)
#}
#
#egc.root <- "/home/gad/Data/EGC/ForGadEGC/genome/egc_1000g_geno05_maf005_hwe6_extract610K_filtered_common_ref"
#system(paste(
#   "plink --noweb", "--bfile", egc.root, "--recodeA",
#   "--extract snps.txt", "--out egc",
#   "--reference-allele ~/Data/DILGOM/Human610/dilgom.ref"
#)) 
#d.egc <- read.table("egc.raw", header=TRUE, stringsAsFactors=FALSE)
#nm.egc <- colnames(d.egc)[-(1:6)]
#x.egc <- as.matrix(d.egc[, -(1:6)])
#y.egc <- egc()
#colnames(x.egc) <- gsub("_[AGCT]", "", colnames(x.egc))
#
## Train on entire DILGOM using selected SNPs
s <- names(rs[rs > 0.7])
#s.common <- s[s %in% colnames(x.egc)]
#b <- drop(ginv(cbind(1, x[, s.common]))) %*% y
#names(b) <- c("(Intercept)", s.common)
#
#p.egc <- cbind(1, x.egc[, s.common]) %*% b
#(r2.egc <- 1 - sum((p.egc - y.egc)^2) / sum((y.egc - mean(y.egc))^2))
#


