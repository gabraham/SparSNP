# Process prediction results

library(ggplot2)
library(glmnet)

# for cross-validation
evalpred.crossval <- function(dir=NULL)
{
   folds <- scan("folds.txt", quiet=TRUE)
   folds <- unique(sort(folds))
   
   if(!exists("uni"))
      uni <- FALSE
   
   cv.auc.d <- pred <- NULL
   cv.auc <- lapply(folds, function(fold) {
      y <- scan(sprintf("y.%02d", fold), quiet=TRUE)
      
      if(length(unique(y)) == 1)
         return(NULL)
   
      if(length(unique(y)) == 1)
         return(NULL)
   
      files <- sort(
         list.files(
   	 pattern=sprintf("^beta\\.csv\\.[[:digit:]]+\\.%02d\\.pred[\\.bz2]*$", fold))
      )
   
      pr <- sapply(files, scan, quiet=TRUE)
   
      lambda <- scan(sprintf("lambda1path.csv.%02d", fold), quiet=TRUE)
      nz <- scan(sprintf("nonzero.csv.%02d", fold), quiet=TRUE)
      
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
         y <- scan(sprintf("multivar_y.%02d", fold), quiet=TRUE)
      
         files <- sort(
            list.files(
   	    pattern=sprintf(
   	       "^multivar_beta.csv.[[:digit:]]+.%02d.pred[\\.bz2]*$", fold))
         )
      
         pr <- sapply(files, scan, quiet=TRUE)
         res <- apply(pr, 2, auc, y=y)
      
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
}

# for testing on validation dataset
evalpred.validation <- function(discovery.dir)
{
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
   
   # glmnet::auc expects y = 0/1
   y <- as.numeric(as.character(factor(y, labels=c(0, 1))))
   res <- apply(pr, 2, auc, y=y)

   d <- data.frame(NonZero=nz, AUC=res)
   cv.auc.d <- d[d$NonZero > 0, ]

   save(cv.auc.d, file="crossval.RData")
}


