# Process prediction results

library(ggplot2)
library(glmnet)
#library(ROCR)

folds <- scan("folds.txt")
folds <- unique(sort(folds))

if(!exists("uni"))
   uni <- FALSE

cv.auc.d <- pred <- NULL
cv.auc <- lapply(folds, function(fold) {
   y <- scan(sprintf("y.%02d", fold))
   
   if(length(unique(y)) == 1)
      return(NULL)

   if(length(unique(y)) == 1)
      return(NULL)

   files <- sort(
      list.files(
	 pattern=sprintf("^beta\\.csv\\.[[:digit:]]+\\.%02d\\.pred[\\.bz2]*$", fold))
   )

   pr <- sapply(files, scan)

   lambda <- scan(sprintf("lambda1path.csv.%02d", fold))
   nz <- scan(sprintf("nonzero.csv.%02d", fold))
   
   # Don't count the intercept, except for first #nonzero
   nz[-1] <- nz[-1] - 1  

   # glmnet::auc expects y = 0/1
   y <- as.numeric(as.character(factor(y, labels=c(0, 1))))
   res <- apply(pr, 2, auc, y=y)

   pred <- list(
      p=lapply(1:ncol(pr), function(j) pr[,j]),
      y=lapply(1:ncol(pr), function(j) y)
   )

   list(
      auc=cbind(
	 lambda=lambda[1:length(res)],
	 AUC=res,
	 NonZero=nz[1:length(res)]
      ),
      prediction=pred
   )
})

pred <- lapply(cv.auc, function(x) x$prediction)
cv.auc.d <- data.frame(do.call(rbind, lapply(cv.auc, function(x) x$auc)))

# remove zeros
cv.auc.d <- cv.auc.d[cv.auc.d$NonZero > 0, ] 

cv.uni.auc.d <- pred.uni <- NULL
if(uni)
{
   cv.uni.auc <- lapply(folds, function(fold) {
      # should be 0/1
      y <- scan(sprintf("multivar_y.%02d", fold))
   
      files <- sort(
         list.files(
	    pattern=sprintf(
	       "^multivar_beta.csv.[[:digit:]]+.%02d.pred[\\.bz2]*$", fold))
      )
   
      pr <- sapply(files, scan)
      res <- apply(pr, 2, auc, y=y)
   
      pred <- list(
         p=lapply(1:ncol(pr), function(j) pr[,j]),
         y=lapply(1:ncol(pr), function(j) y)
      )
   
      # Multivar nonzero doesn't intercept
      nz <- scan(sprintf("multivar_nonzero.csv.%02d", fold))
      nz.thin <- scan(sprintf("multivar_nonzero.csv_thinned.%02d", fold))
      status <- scan(
	 list.files(pattern=sprintf("multivar_[irls|status]+.%02d", fold))
      )

      list(
         auc=cbind(
	    AUC=res,
   	    NonZeroPreThin=nz,
   	    NonZero=nz.thin,
   	    Status=status
         ),
         prediction=pred
      )
   })
   
   pred.uni <- lapply(cv.uni.auc, function(x) x$prediction)
   cv.uni.auc.d <- data.frame(
      do.call(rbind, lapply(cv.uni.auc, function(x) x$auc))
   )
}

save(cv.auc.d, cv.uni.auc.d, pred, pred.uni, file="crossval.RData")

