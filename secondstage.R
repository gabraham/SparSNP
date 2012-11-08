
options(error=dump.frames)

options(stringsAsFactors=FALSE)

library(foreach)
library(methods)
library(doMC)

registerDoMC(cores=10)

ridge <- FMPR::ridge

R2 <- function(pr, y)
{
   1 - sum((pr - y)^2) / sum((y - mean(y))^2)
}

param <- read.table("../params.txt", sep="=", row.names=1,
   stringsAsFactors=FALSE)

root <- param["ROOT", 1]
pheno <- param["PHENO", 1]

nfolds <- 10
folds <- scan("folds.txt")

# reference allele
reff <- "../snps.txt"
ref <- read.table(reff, stringsAsFactors=FALSE)
ref <- cbind(ref, apply(ref, 1, paste, collapse="_"))

dy <- read.table(pheno, stringsAsFactors=FALSE, header=FALSE)

r2 <- foreach(fold= 1:nfolds - 1) %dopar% {

   lf <- list.files(pattern=sprintf("^beta.csv.[[:digit:]]+.%02d$", fold))

   # first one is always an intercept-only model
   lf <- lf[lf != "beta.csv.00.00"]

   topsnps <- lapply(lf, function(f) {
      # read and ignore intercept
      b <- read.table(f, sep=",", row.names=1, header=FALSE)
      as.integer(rownames(b)[-1])
   })

   allsnps <- ref[unique(unlist(topsnps)), 1]
   allsnpsf <- sprintf("allsnps.fold.%d.txt", fold)
   write.table(allsnps, file=allsnpsf, row.names=FALSE,
         col.names=FALSE, quote=FALSE)

   # Extract data (training + test together)
   system(paste(
      "plink --noweb --bfile ", root, " --extract ", allsnpsf,
       " --recodeA --reference-allele ", reff,
      " --out ", allsnpsf, sep=""))

   geno <- read.table(sprintf("%s.raw", allsnpsf), header=TRUE,
      stringsAsFactors=FALSE)
   X <- as.matrix(geno[, -(1:6)])
   if(any(is.na(X)))
      stop("missing values in X")

   if(any(!colnames(X) %in% ref[, 3]))
      stop("SNP name mismatch")

   ytrain <- dy[folds != fold, 3]
   ytest <- dy[folds == fold, 3]

   r <- sapply(topsnps, function(s) {
      cat("n:", length(s), "\n")
      X1train <- cbind(1, X[folds != fold, ref[s, 3]])
      X1test <- cbind(1, X[folds == fold, ref[s, 3]])
      b <- ridge(X1train, ytrain, lambda=1e-3)[[1]]
      pr <- X1test %*% b
      R2(pr, ytest)
   })
   data.frame(Measure=r, NonZero=sapply(topsnps, length))
}

res <- do.call(rbind, r2)

save(res, file="secondstage.RData")

