# Process prediction results

auc <- function(y, p)
{
   r <- rank(p)
   n1 <- sum(y == 1)
   n0 <- length(y) - n1
   u <- sum(r[y == 1]) - n1 * (n1 + 1) / 2
   u / (n0 * n1)
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
	 y <- as.numeric(as.character(factor(y, labels=c(0, 1))))
	 res <- apply(pr, 2, auc, y=y)
      } else {
	 res <- apply(pr, 2, r2, y=y)
      }

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
      
	 pr <- sapply(files, scan, quiet=TRUE)
      
	 if(type == "AUC") {
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
		     fold)), quiet=TRUE
         )

	 w <- which(status == 3)
	 if(length(w) == 0) {
	    w <- length(status)
	 }

	 ok <- which(status[1:w] == 1)
   
         list(
            r=cbind(
	       Measure=res[ok],
	       NonZeroPreThin=nz[ok],
	       NonZero=nz.thin[ok]
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

   ###########################################################################    
   # lasso
   files <- sort(
      list.files(
         pattern="^beta\\.csv\\.[[:digit:]]+\\.[[:digit:]]+\\.pred[\\.bz2]*$")
   )
   
   pr <- sapply(files, scan, quiet=TRUE)
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
   cv.d <- d[d$NonZero > 0, ]

   ###########################################################################    
   # univariable/multivariable
   files <- sort(
      list.files(
         pattern="^multivar_beta\\.csv\\.[[:digit:]]+\\.[[:digit:]]+\\.pred[\\.bz2]*$")
   )
    
   # remove predictions from bed models that didn't converge etc
   lf <- list.files(pattern="multivar_[irls|status]+\\.[[:digit:]]+$",
	 path=discovery.dir, full.names=TRUE)
   status <- lapply(lf, scan, quiet=TRUE)
   names(status) <- lf
   f <- sapply(strsplit(lf, "\\."), tail, n=1)
   l <- lapply(f, function(g) {
      grep(sprintf("^multivar_beta\\.csv\\.[[:digit:]]+\\.%s\\.pred$", g),
	 files, value=TRUE)
   })

   l2 <- lapply(seq(along=l), function(i) {
      l[[i]][status[[i]] == 1]
   })

   l3 <- unlist(l2)

   ## not a sparse model
   pr <- sapply(l3, scan, quiet=TRUE)
   beta <- gsub("\\.pred", "", l3)
   #nz <- sapply(beta, function(x) {
   #   b <- scan(sprintf("%s/%s", discovery.dir, x), quiet=TRUE)
   #   sum(b[-1] != 0)
   #})
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
   cv.uni.d <- d[d$NonZero > 0, ]

   save(cv.d, cv.uni.d, file="crossval.RData")
   list(cv.d, cv.uni.d)
}


