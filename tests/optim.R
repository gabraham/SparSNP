# gradient descent and coordinate descent

softthresh <- function(a, b)
{
   sign(a) * pmax(abs(a) - b, 0)
}

l2loss <- function(x, y, b)
{
   mean((y - x %*% b)^2)
}

l2dloss <- function(x, y, b, j=NULL)
{
   p <- ncol(x)
   if(is.null(j))
      j <- 1:p
   crossprod(x[,j,drop=FALSE], x %*% b - y)
}

l2d2loss <- function(x, b, j=NULL)
{
   p <- ncol(x)
   if(is.null(j))
      j <- 1:p
   diag(crossprod(x[,j]))
}

logloss <- function(x, y, b)
{
   n <- nrow(x)
   s <- sapply(1:n, function(i) {
      log(1 + exp(x[i,] %*% b)) - y[i] * b %*% x[i,]
   })
   mean(s)
}

logdloss <- function(x, y, b, j=NULL)
{
   n <- nrow(x)
   if(is.null(j))
      j <- 1:ncol(x)
   s <- sapply(1:n, function(i) {
      p <- drop(exp(x[i,] %*% b))
      x[i,j,drop=FALSE] * (p / (1 + p) - y[i])
   })
   if(length(j) == 1)
      sum(s)
   else
      rowSums(s)
}

logd2loss <- function(x, b, j=NULL)
{
   n <- nrow(x)
   #if(is.null(j))
   #   j <- 1:ncol(x)
   s <- sapply(1:n, function(i) {
      p <- drop(exp(x[i, ] %*% b))
      P <- p / (1 + p)
      x[i,j,drop=FALSE]^2 * P * (1 - P)
   })
   #if(length(j) == 1)
      sum(s)
   #else
   #   rowSums(s)
}

gd <- function(x, y, lossfunc, dfunc, d2func, maxiter=10)
{
   x <- cbind(1, x)
   p <- ncol(x)
   beta <- numeric(p)
   for(i in 1:maxiter)
   {
      grad <- dfunc(x, y, beta)
      d2 <- d2func(x, beta)
      beta <- beta - grad / d2
      #beta <- beta - 1e-2 * drop(ginv(cov(x)) %*% grad)
      #beta <- beta - crossprod(x, x %*% beta - y) * 1e-6
      cat(i, "loss:", lossfunc(x, y, beta), "\n")
   }
   beta
}

l2d1phi <- function(lp)
{
   lp
}

l2d2phi <- function(lp)
{
   1
}

logd1phi <- function(lp)
{
   1 / (1 + exp(-lp))
}

logd2phi <- function(lp)
{
   p <- exp(lp)
   P <- p / (1 + p)
   P * (1 - P)
}

# naive coordinate descent
cd1 <- function(x, y, lossfunc, dfunc, d2func, maxiter=50)
{
   x <- cbind(1, x)
   n <- nrow(x)
   p <- ncol(x)
   beta <- numeric(p)

   for(i in 1:maxiter)
   {
      for(j in 1:p)
      {
         grad <- dfunc(x, y, beta, j)
         d2 <- d2func(x, beta, j)
	 beta[j] <- beta[j] - grad / d2
      }
      cat(i, "loss:", lossfunc(x, y, beta), "\n")
   }
   beta
}

cd3 <- function(x, y, lossfunc, d1phi, d2phi, lambda1=0, maxiter=20)
{
   x <- cbind(1, x)
   n <- nrow(x)
   p <- ncol(x)
   beta <- numeric(p)
   lp <- x %*% beta

   grad <- d2 <- numeric(p)

   #for(i in 1:maxiter)
   {
      for(j in 1:p)
      {
	 beta.old <- beta[j]

	 grad[j] <- sum(x[, j] * (d1phi(lp) - y))
	 d2[j] <- sum(x[,j]^2 * d2phi(lp))
      }
   }
   list(grad, d2)
}

# coordinate descent
cd2 <- function(x, y, lossfunc, d1phi, d2phi, lambda1=0, maxiter=20)
{
   x <- cbind(1, x)
   n <- nrow(x)
   p <- ncol(x)
   beta <- numeric(p)
   lp <- x %*% beta

   for(i in 1:maxiter)
   {
      cat("epoch", i, "\r")
      for(j in 1:p)
      {
	 beta.old <- beta[j]

	 grad <- sum(x[, j] * (d1phi(lp) - y))
	 d2 <- sum(x[,j]^2 * d2phi(lp))

	 #cat(grad, d2, "\n")

	 if(grad != 0 && d2 != 0)
	 {
	    beta[j] <- if(j > 1) {
	       softthresh(beta.old - grad / d2, lambda1)
	    } else beta.old - grad / d2

	    beta[j] <- sign(beta[j]) * pmin(abs(beta[j]), 20) 
   
	    lp <- lp + x[, j] * (beta[j] - beta.old)
	    lp <- sign(lp) * pmin(abs(lp), 20)
	 }
	 #cat(lp, "\n")
      }
      nz <- sum(beta != 0)
      #cat(i, "loss:", lossfunc(x, y, beta), "nonzero:", nz, "\n")
   }
   cat("\n")
   beta
}


#set.seed(32431431)
#n <- 1000
#p <- 50
#x <- matrix(rnorm(n * p), n, p)
#beta <- rnorm(p + 1)
#y1 <- cbind(1, x) %*% beta + rnorm(n)
#y2 <- ifelse(runif(n) <= plogis(cbind(1, x) %*% beta), 1, 0)
#
#n <- 1e3
#p <- 1e3
#x <- matrix(as.numeric(readBin("sim5/sim.bin", n=n*(p+1), what="raw")),
#      nrow=n, byrow=TRUE)
#y <- x[,1]
#x <- x[,-1]
#
#system.time({b.l2.cd <- cd(x, y, l2loss, l2dloss, l2d2loss)})
#system.time({b.l2.cd2 <- cd2(x, y, l2loss, l2dloss, l2d2loss)})
#system.time({b.l2.gd <- gd(x, y, l2loss, l2dloss, l2d2loss)})
##b.log.gd <- gd(x, y, logloss, logdloss, logd2loss)
##b.log.cd <- cd(x, y, logloss, logdloss, logd2loss)
#
#stop()
#
##b.svd <- ginv(cbind(1, x)) %*% y 
##l2loss(cbind(1, x), y, b.svd)
#
##z <- cbind(y2, x)
##writeBin(as.numeric(t(z)), con="x.bin")
##
#
#f1 <- function(b)
#{
#   l2loss(cbind(1, x), y, b)
#}
#
#g1 <- function(b)
#{
#   l2dloss(cbind(1, x), y, b)
#}
#
#cx <- cov(x)
#
#h1 <- function(b)
#{
#   cx
#}

#system.time({opt1 <- optim(numeric(ncol(x) + 1), f1, g1, method="BFGS")})

#system.time({opt2 <- nlminb(numeric(ncol(x) + 1), f1, g1, h1) })

