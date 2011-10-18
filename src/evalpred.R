# Process prediction results

library(ggplot2)
library(glmnet)

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
   1 - sum((f - y)^2) / sum((y - mean(y))^2)
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
      y <- scan(sprintf("y.%02d", fold), quiet=TRUE)
      
      if(length(unique(y)) == 1)
         return(NULL)
   
      if(length(unique(y)) == 1)
         return(NULL)
   
      files <- sort(
         list.files(
	    pattern=sprintf(
	       "^beta\\.csv\\.[[:digit:]]+\\.%02d\\.pred[\\.bz2]*$", fold))
      )
   
      pr <- sapply(files, scan, quiet=TRUE)
      lambda <- scan(sprintf("lambda1path.csv.%02d", fold), quiet=TRUE)
      nz <- scan(sprintf("nonzero.csv.%02d", fold), quiet=TRUE)
      
      # Don't count the intercept, except for first #nonzero
      nz[-1] <- nz[-1] - 1  
   
      if(type == "AUC") {
	 # glmnet::auc expects y = 0/1
	 y <- as.numeric(as.character(factor(y, labels=c(0, 1))))
	 res <- apply(pr, 2, auc, y=y)
      } else {
	 res <- apply(pr, 2, r2, y=y)
      }

      #browser()

      pred <- list(
         p=lapply(1:ncol(pr), function(j) pr[,j]),
         y=lapply(1:ncol(pr), function(j) y)
      )
   
      list(
         r=cbind(
	    lambda=lambda[1:length(res)],
	    Measure=res,
	    NonZero=nz[1:length(res)]
         ),
         prediction=pred
      )
   })
   
   pred <- lapply(cv, function(x) x$prediction)
   cv.d <- data.frame(do.call(rbind, lapply(cv, function(x) x$r)))
   
   # remove zeros
   cv.d <- cv.d[cv.d$NonZero > 0, ] 
   
   cv.uni.d <- pred.uni <- NULL
   if(uni)
   {
      cv.uni <- lapply(folds, function(fold) {
         # should be 0/1
         y <- scan(sprintf("multivar_y.%02d", fold), quiet=TRUE)
      
         files <- sort(
            list.files(
   	    pattern=sprintf(
   	       "^multivar_beta.csv.[[:digit:]]+.%02d.pred[\\.bz2]*$", fold))
         )
      
	 if(type == "AUC") {
	    # glmnet::auc expects y = 0/1
	    y <- as.numeric(as.character(factor(y, labels=c(0, 1))))
	    res <- apply(pr, 2, auc, y=y)
	 } else {
	    res <- apply(pr, 2, r2, y=y)
	 }

         pred <- list(
            p=lapply(1:ncol(pr), function(j) pr[,j]),
            y=lapply(1:ncol(pr), function(j) y)
         )
      
         # Multivar nonzero doesn't have intercept
         nz <- scan(sprintf("multivar_nonzero.csv.%02d", fold), quiet=TRUE)
         nz.thin <- scan(sprintf("multivar_nonzero.csv_thinned.%02d", fold),
	    quiet=TRUE)
         status <- scan(
	       list.files(pattern=sprintf("multivar_[irls|status]+.%02d",
		     fold), quiet=TRUE)
         )
   
         list(
            r=cbind(
	       Measure=res,
	       NonZeroPreThin=nz,
	       NonZero=nz.thin,
	       Status=status
            ),
            prediction=pred
         )
      })
      
      pred.uni <- lapply(cv.uni, function(x) x$prediction)
      cv.uni.d <- data.frame(
         do.call(rbind, lapply(cv.uni, function(x) x$r))
      )
   }
   
   save(cv.d, cv.uni.d, pred, pred.uni, file="crossval.RData")
}

# for testing on validation dataset
evalpred.validation <- function(type=NULL, discovery.dir=NULL)
{
   type <- match.arg(type, c("AUC", "R2"))
   if(is.null(discovery.dir))
      stop("discovery.dir not specified")

   y <- scan("y.txt", quiet=TRUE)

   files <- sort(
      list.files(
         pattern="^beta\\.csv\\.[[:digit:]]+\\.[[:digit:]]+\\.pred[\\.bz2]*$")
   )
   
   pr <- sapply(files, scan, quiet=TRUE)

   cat("pwd:", getwd(), "\n")

   beta <- gsub("\\.pred", "", files)
   nz <- sapply(beta, function(x) {
      nrow(read.table(sprintf("%s/%s", discovery.dir, x))) - 1
   })
   
   if(type == "AUC") {
      # glmnet::auc expects y = 0/1
      y <- as.numeric(as.character(factor(y, labels=c(0, 1))))
      res <- apply(pr, 2, auc, y=y)
   } else {
      res <- apply(pr, 2, r2, y=y)
   }

   d <- data.frame(NonZero=nz, Measure=res)
   cv.d <- d[d$NonZero > 0, ]

   save(cv.d, file="crossval.RData")
}


