library(glmnet)

options(digits=9)

#set.seed(439)

#n <- 100
#p <- 1000
#x <- matrix(sample(0:2, n * p, replace=TRUE), n, p)
#
## scale to unit norm, not unit variance
#x <- apply(x, 2, function(z) {
#   z2 <- z - mean(z) 
#   z2 / sqrt(mean(z2^2))
#})
#
#x1 <- cbind(1, x)
#w <- sample(0:1, p + 1, replace=TRUE, prob=c(0.9, 0.1))
#wl <- as.logical(w)
#beta <- rnorm(p + 1, 0, 0.5) * w
#y <- drop(x1 %*% beta) + rnorm(n)

soft <- function(b, g)
{
   sign(b) * pmax(abs(b) - g, 0)
}

cd.plain <- function(x, y, lambda1=0, maxepoch=100, maxiter=50,
      eps=1e-9)
{
   p <- ncol(x)
   n <- nrow(x)
   lp <- numeric(n)
   beta <- numeric(p)
   
   for(epoch in 1:maxepoch)
   {
      for(j in 1:p)
      {
	 for(iter in 1:maxiter)
	 {
	    s <- drop(crossprod(x[,j], lp - y)) / n
	    beta_new <- beta[j] - s
	    if(j > 1)
	       beta_new <- soft(beta_new, lambda1)

	    diff <- beta_new - beta[j]
	    lp <- lp + x[, j] * diff
	    beta[j] <- beta_new
	    if(abs(diff) < eps)
	       break
	 }
      }
   }
   beta
}

cd.activeset <- function(x, y, lambda1=0, beta=numeric(ncol(x)),
   lp=numeric(nrow(x)), maxepoch=500, maxiter=50, eps=1e-4)
{
   p <- ncol(x)
   n <- nrow(x)
   beta_old <- beta
   active_new <-  rep(TRUE, p)
   active_old <- rep(FALSE, p)
   allconverged <- 0

   cat("lambda1:", lambda1, "\n")
   for(epoch in 1:maxepoch)
   {
      for(j in 1:p)
      {
	 if(active_new[j]) 
	 {
	    for(iter in 1:maxiter)
	    {
	       s <- sum(x[,j] * (lp - y)) / n
               beta_new <- beta[j] - s
               if(j > 1)
                  beta_new <- soft(beta_new, lambda1)

	       diff <- beta_new - beta[j]
	       lp <- lp + x[, j] * diff
	       beta[j] <- beta_new
	       if(abs(diff) < eps)
	          break
	    }
	    if(iter == maxiter)
	       cat("maxiter reached for variable", j, "\n")
	    active_new[j] <- beta[j] != 0
	 }
      }

      if(all(abs(beta - beta_old) <= eps))
      {
         allconverged <- allconverged + 1

         if(allconverged == 1) {
            active_old <- active_new
            active_new[] <- TRUE
         } else if(allconverged == 2) {
            if(all(active_old == active_new))
               break
            active_old <- active_new
            active_new[] <- TRUE
            allconverged <- 1
         }
      }
      else
	 allconverged <- 0 

      beta_old <- beta
   }

   if(epoch == maxepoch)
      cat("failed to converged within", maxepoch, "epochs\n")

   cat("\n")
   beta
}

# Slightly different formulation
cd.activeset2 <- function(x, y, lambda1=0, lambda1max=0, beta=numeric(ncol(x)),
   lp=numeric(nrow(x)), maxepoch=500, maxiter=50, eps=1e-4)
{
   p <- ncol(x)
   n <- nrow(x)
   beta_old <- beta
   active_new <- rep(TRUE, p)
   active_old <- rep(FALSE, p)
   allconverged <- 0
   r <- lp - y
   x2 <- abs(drop(crossprod(x, r)))
   never <- x2 < 2 * lambda1 - lambda1max
   active_new[never] <- FALSE

   cat("never:", sum(never), "\n")
   
   cat("lambda1:", lambda1, "\n")
   for(epoch in 1:maxepoch)
   {
      for(j in 1:p)
      {
	 if(active_new[j])
	 {
	    for(iter in 1:maxiter)
	    {
	       s <- sum(x[,j] * (lp - y)) / n
	       beta_new <- beta[j] - s
	       if(j > 1)
		  beta_new <- soft(beta_new, lambda1)

	       diff <- beta_new - beta[j]
	       lp <- lp + x[, j] * diff
	       beta[j] <- beta_new
	       r <- lp - y

	       if(abs(diff) < eps)
	          break
	    }
	    if(iter == maxiter)
	       cat("maxiter reached for variable", j, "\n")
	    active_new[j] <- beta[j] != 0

	    x2 <- abs(drop(crossprod(x, r)))
	    never <- x2 < 2 * lambda1 - lambda1max
	    active_new[never] <- FALSE
	 }
      }

      if(all(abs(beta - beta_old) <= eps))
      {
         allconverged <- allconverged + 1

         if(allconverged == 1) {
            active_old <- active_new
            active_new[!never] <- TRUE
         } else if(allconverged == 2) {
            if(all(active_old == active_new))
               break
            active_old <- active_new
            active_new[!never] <- TRUE
            allconverged <- 1
         }
      }
      else
	 allconverged <- 0 

      beta_old <- beta
   }

   if(epoch == maxepoch)
      cat("failed to converged within", maxepoch, "epochs\n")

   cat("\n")
   beta
}


#system.time({
#   g <- glmnet(x, y, nlambda=20)
#})
#b.glmnet <- as.matrix(coef(g))
#nlambda <- length(g$lambda)
#
#matplot(g$lambda, t(b.glmnet), type="l", lty=1, col=1, lwd=5,
#      log="x", pch=21)
#points(rep(min(g$lambda), p+1), coef(lm(y ~ x)))
#
#b.cd.as <- matrix(0, p + 1, nlambda + 1)
#system.time({
#   for(i in 1:nlambda + 1)
#   {
#      b <- b.cd.as[, i - 1]
#      lp <- drop(x1 %*% b)
#      b.cd.as[, i] <- cd.activeset(x=x1, y=y, beta=b,
#	    lp=lp, lambda1=g$lambda[i-1])
#   }
#})
#b.cd.as <- b.cd.as[, -1]
#
#matlines(g$lambda, t(b.cd.as), type="l", lty=1, col=2, pch=21, lwd=3)
#
#b.cd.as2 <- matrix(0, p + 1, nlambda + 1)
#system.time({
#   for(i in 1:nlambda + 1)
#   {
#      b <- b.cd.as2[, i - 1]
#      lp <- drop(x1 %*% b)
#      b.cd.as2[, i] <- cd.activeset2(x=x1, y=y, beta=b,
#	    lp=lp, lambda1=g$lambda[i-1], lambda1max=g$lambda[1])
#   }
#})
#b.cd.as2 <- b.cd.as2[, -1]
#
#matlines(g$lambda, t(b.cd.as2), type="l", lty=1, col=3, pch=21, lwd=3)
#sqrt(mean((b.glmnet - b.cd.as2)^2))

#system.time({
#   b.cd.plain <- sapply(g$lambda, cd.plain, x=x1, y=y)
#})
#matlines(g$lambda, t(b.cd.plain), type="l", lty=1, col=7, pch=21)
#sqrt(mean((b.glmnet - b.cd.plain)^2))

