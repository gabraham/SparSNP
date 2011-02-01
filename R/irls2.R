
library(glmnet)
library(MASS)

irls <- function(x, y, lambda=0, dispersion=1, maxiter=250, scale=FALSE,
      maxlp=700, thresh=1e-3)
{
   if(NROW(x) != length(y))
      stop("NROW(x) != length(y)")

   x <- cbind(x)
   good <- apply(x, 1, function(x) all(!is.na(x)))
   x <- x[good,]
   x1 <- if(scale) {
      cbind(1, scale(x))
   } else {
      cbind(1, x)
   }
   y <- y[good]

   if(sum(good) == 0)
      stop("all samples have missing observations")

   cat("good:", sum(good), "\n")
   
   converged <- FALSE
   n <- length(y)
   p <- ncol(x1)
   lp <- numeric(n)
   lpinv <- 1 / (1 + exp(-lp))
   beta <- numeric(p)
   #s <- 10
   iter <- 1
   #L <- diag(c(0, rep(lambda, p - 1)))

   cat(dim(x), "\n")
   dev <- Inf

   while(iter <= maxiter)
   {
      cat("iter", iter, "dev:", dev, "\n")
      lp <- drop(x1 %*% beta)
      #lp <- sign(lp) * pmin(abs(lp), maxlp)
      lpinv <- 1 / (1 + exp(-lp))
      w <- lpinv * (1 - lpinv)

      grad <- crossprod(x1, lpinv - y)
      hess <- crossprod(x1, diag(w)) %*% x1# + L
      invhess <- ginv(hess)
      beta <- beta - invhess %*% grad

      #z <- lp + (1 / w) * (y - lpinv)
      #l <- try(lm(z ~ x, weights=w))
      #if(class(l) == "try-error")
      #   browser()
      #betaold <- beta
      #beta <- coef(l)
      #beta[is.na(beta)] <- 0
      #s <- beta - betaold

      if(iter > 1) {
	 dev.old <- dev
      }

      dev <- mean(log(1 + exp(lp)) - y * lp) #+ lambda * sum(beta^2)

      browser()
      if(any(abs(beta) > 30) || dev <= 1e-2)
      {
	 converged <- FALSE
	 warning("irls diverged\n")
	 break
      }
      if(iter > 1 && abs(dev - dev.old) / (abs(dev) + 0.1) < thresh)
      {
	 converged <- TRUE
	 break
      }

      iter <- iter + 1
   }


   hess <- crossprod(x1, diag(w)) %*% x1
   invhess <- ginv(hess)
   #invhess <- solve(hess)
   se <- sqrt(diag(invhess) * dispersion)
   zscore <- beta / se
   pval <- 2 * pnorm(abs(zscore), lower.tail=FALSE)
   
   res <- cbind(
      coef=as.numeric(beta),
      se=se,
      zscore=as.numeric(zscore),
      pvalue=as.numeric(pval))

   list(coef=res, dev=dev, hessian=hess, lambda=lambda, converged=converged)
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

