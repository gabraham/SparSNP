

l2loss <- function(x, y, b)
{
   sum((y - x %*% b)^2)
}

l2dloss <- function(x, y, b, j=NULL)
{
   p <- ncol(x)
   if(is.null(j))
      j <- 1:p
   crossprod(x[,j,drop=FALSE], x %*% b - y)
}

l2d2loss <- function(x, y, b, j=NULL)
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
   sum(s)
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

logd2loss <- function(x, y, b, j=NULL)
{
   n <- nrow(x)
   if(is.null(j))
      j <- 1:ncol(x)
   s <- sapply(1:n, function(i) {
      p <- drop(exp(x[i, ] %*% b))
      x[i,j,drop=FALSE]^2 * p^2 / (1 + p)^2
   })
   if(length(j) == 1)
      sum(s)
   else
      rowSums(s)
}

gd <- function(x, y, lossfunc, dfunc, d2func, maxiter=10)
{
   x <- cbind(1, x)
   p <- ncol(x)
   beta <- numeric(p)
   for(i in 1:maxiter)
   {
      grad <- dfunc(x, y, beta)
      d2 <- d2func(x, y, beta)
      beta <- beta - grad / d2
      cat(i, "loss:", lossfunc(x, y, beta), "\n")
   }
   beta
}

cd <- function(x, y, lossfunc, dfunc, d2func, maxiter=10)
{
   x <- cbind(1, x)
   p <- ncol(x)
   beta <- numeric(p)
   for(i in 1:maxiter)
   {
      for(j in 1:p)
      {
         grad <- dfunc(x, y, beta, j)
         d2 <- d2func(x, y, beta, j)
         beta[j] <- beta[j] - grad / d2
      }
      cat(i, "loss:", lossfunc(x, y, beta), "\n")
   }
   beta
}

set.seed(32431431)
n <- 1000
p <- 50
x <- matrix(rnorm(n * p), n, p)
beta <- rnorm(p + 1)
y1 <- cbind(1, x) %*% beta + rnorm(n)
y2 <- ifelse(runif(n) <= plogis(cbind(1, x) %*% beta), 1, 0)

b.l2.gd <- gd(x, y1, l2loss, l2dloss, l2d2loss)
b.l2.cd <- cd(x, y1, l2loss, l2dloss, l2d2loss)
b.log.gd <- gd(x, y2, logloss, logdloss, logd2loss)
b.log.cd <- cd(x, y2, logloss, logdloss, logd2loss)

