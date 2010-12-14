
library(glmnet)
library(MASS)

irls <- function(x, y, lambda=0, dispersion=1, maxiter=250)
{
   good <- apply(x, 1, function(x) all(!is.na(x)))
   x <- x[good,]
   x1 <- cbind(1, x)
   y <- y[good]

   if(sum(good) == 0)
      stop("all samples have missing observations")

   cat("good:", sum(good), "\n")
   
   n <- length(y)
   p <- ncol(x1)
   lp <- numeric(n)
   lpinv <- 1 / (1 + exp(-lp))
   beta <- numeric(p)
   s <- 10
   iter <- 1
   L <- diag(c(0, rep(lambda, p - 1)))


   cat(dim(x), "\n")

   while(any(abs(s) >= 1e-9) && iter <= maxiter)
   {
      cat("iter", iter, "\n")
      lp <- drop(x1 %*% beta)
      lpinv <- 1 / (1 + exp(-lp))
      w <- lpinv * (1 - lpinv)

      #grad <- crossprod(x1, lpinv - y)
      #cat("grad:", grad, "\n")
      #hess <- crossprod(x1, diag(w)) %*% x1 + L
      #cat("hess:", hess, "\n")
      #invhess <- ginv(hess)
      #cat("invhess:", invhess, "\n")
      #s <- invhess %*% grad
      #beta <- beta - s

      z <- lp + (1 / w) * (y - lpinv)
      l <- try(lm(z ~ x, weights=w))
      if(class(l) == "try-error")
	 return(NULL)

      betaold <- beta
      beta <- coef(l)
      bad <- is.na(beta)
      s <- beta[!bad] - betaold[!bad]
      iter <- iter + 1
   }

   hess <- crossprod(x1, diag(w)) %*% x1
   invhess <- ginv(hess)
   se <- sqrt(diag(invhess) * dispersion)
   zscore <- beta / se
   pval <- 2 * pnorm(abs(zscore), lower.tail=FALSE)
   
   res <- cbind(
      coef=as.numeric(beta),
      se=se,
      zscore=as.numeric(zscore),
      pvalue=as.numeric(pval))
   res
}

cv.irls <- function(x, y, nfolds=10, nreps=10)
{
   n <- length(y)
   sapply(1:nreps, function(i) {
      folds <- sample(1:nfolds, n, TRUE)
      sapply(1:nfolds, function(fold) {
         res <- irls(x[folds != fold, ], y[folds != fold], maxiter=20)
         pr <- cbind(1, x[folds == fold, ]) %*% res[,1]
         auc(y[folds == fold], pr)
      })
   })
}

