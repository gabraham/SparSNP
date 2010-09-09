library(glmnet)

options(digits=9)

set.seed(123456)

n <- 100
p <- 30
x <- scale(matrix(rnorm(n * p), n, p))
x1 <- cbind(1, x)
beta <- rnorm(p + 1, 0, 0.5)
y <- drop(x1 %*% beta) + rnorm(n)

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
   lp=numeric(nrow(x)), maxepoch=100, maxiter=50, eps=1e-4, ztol=1e-9)
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
	       s <- drop(crossprod(x[,j], lp - y)) / n
               beta_new <- beta[j] - s
               if(j > 1)
                  beta_new <- soft(beta_new, lambda1)

	       if(abs(beta_new) <= ztol)
		  beta_new <- 0

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

      beta_old <- beta
   }

   if(epoch == maxepoch)
      cat("failed to converged within", maxepoch, "epochs\n")

   cat("\n")
   beta
}

g <- glmnet(x, y, nlambda=20)
b.glmnet <- as.matrix(coef(g))
nlambda <- length(g$lambda)

matplot(g$lambda, t(b.glmnet), type="l", lty=1, col=1, lwd=5,
      log="x", pch=21)
points(rep(min(g$lambda), p+1), coef(lm(y ~ x)))

b.cd.as <- matrix(0, p + 1, nlambda + 1)
system.time({
   for(i in 1:nlambda + 1)
   {
      b <- b.cd.as[, i - 1]
      lp <- drop(x1 %*% b)
      b.cd.as[, i] <- cd.activeset(x=x1, y=y, beta=b,
	    lp=lp, lambda1=g$lambda[i-1])
   }
})
b.cd.as <- b.cd.as[, -1]

matlines(g$lambda, t(b.cd.as), type="l", lty=1, col=4, pch=21, lwd=3)
sqrt(mean((b.glmnet - b.cd.as)^2))

system.time({
   b.cd.plain <- sapply(g$lambda, cd.plain, x=x1, y=y)
})
matlines(g$lambda, t(b.cd.plain), type="l", lty=1, col=7, pch=21)
sqrt(mean((b.glmnet - b.cd.plain)^2))

