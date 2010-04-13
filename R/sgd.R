
#setGeneric("train", function(object, x, y, ...) standardGeneric("train"))
#setGeneric("predict")
#setGeneric("crossval", function(object, x, y, ...) standardGeneric("crossval"))
#setGeneric("bootstrap", function(object, x, y, ...) standardGeneric("bootstrap"))
#setGeneric("bag", function(object, x, y, ...) standardGeneric("bag"))
#
#setClass("sgd", representation(
#  B="matrix", model="character", lambda1="numeric",
#  lambda2="numeric", features="integer", subset="integer",
#  p="numeric", epochs="integer", nsamples="integer", source="character",
#  losses="numeric", stepsize="numeric", anneal="numeric")
#), prototype())

################################################################################
# Loss functions
# 
# These functions return the *gradient* wrt b, not the derivative.

# Squared loss for regression
l2loss <- function(x, y, b)
{
   sum((y - x %*% b)^2)
}

l2dloss <- function(x, y, b)
{
   #crossprod(y - x %*% b, -x)
   crossprod(x, x %*% b - y)
}

# Hinge loss for SVM classification
hingeloss <- function(x, y, b)
{
   sum(pmax(0, 1 - y * x %*% b))
}

hingedloss <- function(x, y, b)
{
   I <- as.numeric(drop(1 - y * x %*% b > 0))
   crossprod(-I * x, y)
}

# Binary logistic, negative of the binomial log likelihood
logloss <- function(x, y, b)
{
   sum(-y * x %*% b + log(1 + exp(x %*% b)))
}

# Note that binary logistic and multinomial logistic are parameterised here
# differently, so a 2-class multinomial will not give you the same
# estimated coefficients as a binary logistic!
logdloss <- function(x, y, b)
{
   p <- exp(x %*% b) 
   crossprod(x, p / (1 + p) - y)
}

# Multinomial logistic
# y is an indicator matrix
# b is a matrix of coefficients
mlogloss <- function(x, Y, B)
{
   p <- x %*% B
   -sum(Y * (p - log(rowSums(exp(p)))))
}

mlogdloss <- function(x, Y, B)
{
   q <- exp(x %*% B)
   p <- q / rowSums(q)
   crossprod(x, p - Y)
}

# Cox proportional hazards
coxloss <- function(x, y, b)
{
   
}

coxdloss <- function(x, y, b)
{
}


#################################################################################
#
#test.multinom <- function()
#{
#   set.seed(10278)
#   K <- 2
#   N <- 1000
#   p <- 10
#   x <- matrix(rnorm(N * p), N, p)
#   beta <- rnorm(p)
#   #y <- sample(1:K, N, TRUE)
#   y <- as.numeric(factor(sign(x %*% beta)))
#   #g <- glmnet(x, factor(y), family="multinomial", alpha=0, lambda=0)
#   #pr <- predict(g, x)
#   Y <- sapply(1:K, function(i) as.numeric(y == i))
#   stepsize <- 1e-3
#
#   maxepochs <- 100
#   B <- matrix(0, p, K)
#   loss <- numeric(N)
#   for(epoch in 1:maxepochs)
#   {
#      for(i in 1:N)
#      {
#	 loss[i] <- mlogloss(x[i, , drop=FALSE], Y[i,,drop=FALSE], B) 
#	 #cat("loss:", loss[i], "\n")
#	 B <- B - stepsize * mlogdloss(x[i,,drop=FALSE], Y[i,,drop=FALSE], B)
#      }
#   }
#
#   q <- exp(x %*% B)
#   pr1 <- q / rowSums(q)
#   mean(apply(pr1, 1, which.max) != y)
#   #mean(apply(pr2, 1, which.max) != y)
#}

################################################################################
# Generics and definitions

sgd <- function(x, y, ...)
{
   UseMethod("sgd")
}

rfesgd <- function(x, y, p, nfeats, ...)
{
   UseMethod("rfesgd")
}

setClass("gd",
   representation(B="matrix", model="character",
      lambda1="numeric", lambda2="numeric", lambdaE="numeric", alpha="numeric",
      threshold="numeric", stepsize="numeric", anneal="numeric", subset="logical",
      features="integer", loss="numeric"),
   prototype(B=matrix(), model=character(), lambda1=numeric(),
      lambda2=numeric(), lambdaE=numeric(), alpha=numeric(),
      threshold=numeric(), stepsize=numeric(), anneal=numeric(), subset=logical(),
      features=integer(), loss=numeric())
)

setClass("sgd", contains="gd")

setClass("sgdMem", contains="sgd")

setClass("sgdDisk", contains="sgd")

setMethod("coef", signature(object="sgd"), function(object) object@B)
setMethod("coefficients", signature(object="sgd"), function(object) object@B)

################################################################################
# Fitting functions
#

sgd.gmatrix <- function(g, B=NULL,
      model=c("linear", "logistic", "hinge", "multinomial"),
      lambda1=0, lambda2=0, lambdaE=0, alpha=numeric(), threshold=1e-3,
      stepsize=1e-3, maxepochs=50, anneal=stepsize, blocksize=1,
      maxiter=Inf, subset=logical(), saveloss=FALSE, scale=list(mean=0, sd=1),
      features=1:g@ncol,
      verbose=TRUE, ylevels=NULL)
{
   if(length(subset) > 1 && blocksize > 1)
      stop("blocksize > 1 currently not supported when subset is used")

   if(length(subset) == 0)
      subset <- rep(TRUE, g@nrow)

   model <- match.arg(model)
   loss <- switch(model,
	 linear=l2loss,
	 logistic=logloss,
	 hinge=hingeloss,
	 multinomial=mlogloss
   )
   dloss <- switch(model,
	 linear=l2dloss,
      	 logistic=logdloss,
      	 hinge=hingedloss,
      	 multinomial=mlogdloss
   )

   nf <- length(features)
   scale$mean <- if(length(scale$mean) == g@ncol) {
      scale$mean[features]
   } else scale$mean

   scale$sd <- if(length(scale$sd) == g@ncol) {
      scale$sd[features] 
   } else scale$sd

   # Don't penalise intercept
   lambda1 <- c(0, rep(lambda1, nf))
   lambda2 <- c(0, rep(lambda2, nf))

   if(model == "multinomial") {
      classes <- sort(unique(drop(g@companions$y)))
      K <- length(classes)
      if(is.null(B)) {
	 B <- matrix(0, nf + 1, K)
	 colnames(B) <- 1:ncol(B)
      }
      g@companions$y <- sapply(classes, function(i) {
	 as.numeric(g@companions$y == i)
      })
   } else {
      if(is.null(B))
	 B <- matrix(0, nf + 1, 1)
      #names(b) <- names(b.best) <- 1:length(b)
   }
   B.best <- B

   losses <- rep(0, maxepochs + 1)
   epoch <- epoch.best <- 2

   # Loop over epochs
   while(TRUE)
   {
      # Loop over samples
      for(i in 1:g@nrow)
      {
	 r <- nextRow(g)
	 if(length(r) == 0 || i > maxiter)
	    break

	 if(subset[i]) {
	    x <- r[[1]]
	    y <- r[[2]]$y
	    x <- cbind(1, (x - scale$mean) / scale$sd)
	    l <- loss(x, y, B)
	    
	    ## Ignore samples that make NaN loss (especially relevant for log
	    ## loss) 
	    #if(is.nan(l)) {
	    #   stepsize <- stepsize / 2
	    #} else {
	    cat(i, "sample loss:", l, "\r")
	       losses[epoch] <- losses[epoch] + l
	       grad <- dloss(x, y, B) + lambda2 * B + lambda1 * sign(B)
	       B <- B - stepsize * grad
	    #}
	 } else if(verbose) {
	    cat("skipping", i, "\n")
	 }
      }
      cat("\n")

      # Step halving and greedy choice of best parameters
      if(epoch > 2 && losses[epoch] > losses[epoch-1]) {
         stepsize <- stepsize / 2 
	 B <- B.best
	 if(verbose)
	    cat("Reduced step size\n")
      } else {
         stepsize <- stepsize / (1 + anneal)
         B.best <- B
         epoch.best <- epoch
      }
         
      if(verbose) {
         cat("Epoch", epoch-1, ", loss:", losses[epoch], "diff:",
               losses[epoch-1] - losses[epoch], 
               "stepsize:", stepsize,
               "\n")
      }

      # Check for convergence
      if(epoch >= maxepochs
	 || abs(losses[epoch-1] - losses[epoch]) < threshold) {
	 break
      }
      epoch <- epoch + 1
   }

   if(verbose) {
      cat("Best solution at Epoch", epoch.best - 1,
	 "loss:", losses[epoch.best], "\n")
   }

   new("sgd", B=B.best, model=model, lambda1=lambda1, lambda2=lambda2,
	 lambdaE=lambdaE, alpha=alpha, subset=subset, features=features,
	 stepsize=stepsize, anneal=anneal,
	 loss=if(saveloss) losses[-1][2:epoch-1] else as.numeric(NA)
   )
}

# Batch gradient descent
gd <- function(x, y, model=c("linear", "logistic", "hinge", "multinomial"),
      lambda1=0, lambda2=0, lambdaE=0, alpha=NULL,
      threshold=1e-3, stepsize=1 / nrow(x), anneal=stepsize, maxiter=100,
      scale=TRUE)
{
   model <- match.arg(model)
   loss <- switch(model, linear=l2loss, logistic=logloss, hinge=hingeloss,
	 multinomial=mlogloss)
   dloss <- switch(model, linear=l2dloss, logistic=logdloss, hinge=hingedloss,
	 multinomial=mlogdloss)

   K <- length(unique(y))
   x <- if(scale) {
      cbind(1, scale(x))
   } else cbind(1, x)
   p <- ncol(x)
   if(model == "multinomial") {
      b <- matrix(0, p, K)
      y <- sapply(1:K, function(i) as.numeric(y == i))
   } else {
      b <- matrix(0, p, 1)
   }
   l.old <- Inf
   l.new <- 0

   # Elastic net
   if(lambdaE > 0 && length(alpha) > 0) {
      lambda1 <- lambdaE * alpha
      lambda2 <- lambdaE * (1 - alpha)
   }

   # Don't penalise intercept
   lambda1 <- c(0, rep(lambda1, p - 1))
   lambda2 <- c(0, rep(lambda2, p - 1))

   i <- 1
   while(i <= maxiter && max(abs(l.old - l.new)) >= threshold)
   {
      l.old <- l.new
      grad <- dloss(x, y, b) + lambda2 * b + lambda1 * sign(b)
      b <- b - stepsize * grad
      l.new <- loss(x, y, b)
      cat(l.new, "\n")
      stepsize <- stepsize / (1 + anneal)
      i <- i + 1
   }
   cat("converged at", i-1, "iterations\n")

   structure(b, class="gd", model=model, iter=i)
}

#sgd.matrix <- function(x, y,
#      model=c("linear", "logistic", "hinge", "multinomial"),
#      lambda1=0, lambda2=0, lambdaE=0,
#      alpha=NULL, threshold=1e-3, stepsize=1 / nrow(x),
#      maxepochs=1, anneal=stepsize, maxiter=Inf, scale=TRUE)
#{
#   model <- match.arg(model)
#   loss <- switch(model, linear=l2loss, logistic=logloss, hinge=hingeloss,
#	 multinomial=mlogloss)
#   dloss <- switch(model, linear=l2dloss, logistic=logdloss, hinge=hingedloss,
#	 multinomial=mlogdloss)
#
#   x <- if(scale) {
#      cbind(1, scale(x))
#   } else {
#      cbind(1, x)
#   }
#   p <- ncol(x)
#   n <- nrow(x)
#   n <- min(maxiter, n)
#
#   # Don't penalise intercept
#   lambda1 <- c(0, rep(lambda1, p - 1))
#   lambda2 <- c(0, rep(lambda2, p - 1))
#   
#   if(model == "multinomial") {
#      K <- length(unique(y))
#      b <- b.best <- matrix(0, nf + 1, K)
#      colnames(b) <- colnames(b.best) <- 1:ncol(b)
#      y <- sapply(1:K, function(i) as.numeric(y == i))
#   } else {
#      b <- b.best <- rep(0, nf + 1)
#      names(b) <- names(b.best) <- 1:length(b)
#   }
#
#   l.old <- Inf
#   l.new <- 0
#   epoch <- 1
#   while(epoch <= maxepochs && max(abs(l.old - l.new)) >= threshold)
#   {
#      l.new <- l.old <- 0
#      for(i in 1:n)
#      {
#	 grad <- (dloss(x[i, , drop=FALSE], y[i], b)
#	       + lambda2 * b + lambda1 * sign(b))
#	 l <- loss(x[i, , drop=FALSE], y[i], b) 
#	 l.new <- l.new + l
#	 b <- b - stepsize * grad
#      }
#      stepsize <- stepsize / (1 + anneal)
#      cat("Epoch", epoch, ", loss:", l.new, "\n")
#      epoch <- epoch + 1
#   }
#   structure(b, class="sgd", model=model, p=p, epochs=epoch, source="matrix") 
#}

#sgd.character <- function(x="", y, p,
#      model=c("linear", "logistic", "hinge", "multinomial"),
#      lambda1=0, lambda2=0, lambdaE=0, alpha=NULL, threshold=1e-3,
#      stepsize=1e-4, maxepochs=50, anneal=stepsize, blocksize=1,
#      maxiter=Inf, subset=NULL, saveloss=FALSE, mu=0, norm=1, features=1:p,
#      verbose=TRUE)
#{
#   if(nchar(x) == 0)
#      stop("filename x not supplied")
#
#   if(!is.null(subset) && blocksize > 1)
#      stop("blocksize > 1 currently not supported when subset is used")
#
#   if(is.null(subset))
#      subset <- rep(TRUE, length(y))
#
#   model <- match.arg(model)
#   loss <- switch(model,
#	 linear=l2loss,
#	 logistic=logloss,
#	 hinge=hingeloss,
#	 multinomial=mlogloss
#   )
#   dloss <- switch(model,
#	 linear=l2dloss,
#      	 logistic=logdloss,
#      	 hinge=hingedloss,
#      	 multinomial=mlogdloss
#   )
#
#   fname <- x
#   f <- file(fname, open="rb")
#
#   nf <- length(features)
#   if(length(mu) == p)
#      mu <- mu[features]
#
#   if(length(norm) == p)
#      norm <- norm[features] 
#
#   # Don't penalise intercept
#   lambda1 <- c(0, rep(lambda1, nf))
#   lambda2 <- c(0, rep(lambda2, nf))
#
#   if(model == "multinomial") {
#      K <- length(unique(y))
#      b <- b.best <- matrix(0, nf + 1, K)
#      colnames(b) <- colnames(b.best) <- 1:ncol(b)
#      y <- sapply(1:K, function(i) as.numeric(y == i))
#   } else {
#      b <- b.best <- rep(0, nf + 1)
#      names(b) <- names(b.best) <- 1:length(b)
#   }
#
#   #b <- b.best <- rep(0, nf + 1)
#   #names(b) <- names(b.best) <- 1:length(b)
#   losses <- rep(0, maxepochs + 1)
#   epoch <- epoch.best <- 2
#   #cat("features:", features, "\n")
#
#   # Loop over epochs
#   while(TRUE)
#   {
#      # Loop over samples
#      i <- 1
#      while(TRUE)
#      {
#	 if(verbose)
#	    cat(i, "reading data... ")
#	 dat <- readBin(f, what="numeric", n=p * blocksize)
#	 if(length(dat) == 0 || i > maxiter)
#	    break
#	 if(verbose)
#	    cat("read", length(dat), "items; ")
#
#	 if(subset[i]) {
#	    x <- matrix(dat, nrow=min(blocksize, length(dat) / p), ncol=p,
#	          byrow=TRUE)[, features, drop=FALSE]
#	    x <- cbind(1, (x - mu) / norm)
#	    k <- ((i-1) * blocksize + 1):(i * blocksize)
#	    k <- k[1:nrow(x)]
#	    yk <- if(model == "multinomial") {
#	       y[k, , drop=FALSE]
#	    } else y[k]
#	    l <- loss(x, yk, b)
#	    l <- if(is.nan(l)) {
#	       stepsize <- stepsize / 2
#	       0
#	    } else{
#	       l
#	    }
#	    losses[epoch] <- losses[epoch] + l
#	    if(verbose)
#	       cat("loss:", losses[epoch], "\n")
#	    grad <- (dloss(x[i, , drop=FALSE], y[i], b)
#		  + lambda2 * b + lambda1 * sign(b))
#	    if(verbose > 1)
#	       cat(i, grad[1:10], "\n") 
#	    b <- b - stepsize * grad
#	    
#	 } else if(verbose) {
#	    cat("skipping", i, "\n")
#	 }
#	 i <- i + 1
#      }
#
#      if(losses[epoch] > losses[epoch-1]) {
#	 stepsize <- stepsize / 2 
#      } else {
#	 stepsize <- stepsize / (1 + anneal)
#	 b.best <- b
#	 epoch.best <- epoch
#      }
#	 
#      if(verbose) {
#	 cat("Epoch", epoch-1, ", loss:", losses[epoch], "diff:",
#	       abs(losses[epoch-1] - losses[epoch]), 
#	       "stepsize:", stepsize,
#	       "\n")
#      }
#
#      if(epoch == maxepochs
#	 || abs(losses[epoch-1] - losses[epoch]) < threshold) {
#	 #cat("converged", abs(losses[epoch-1] - losses[epoch]) < threshold, "\n")
#	 break
#      }
#      epoch <- epoch + 1
#      close(f)
#      f <- file(fname, open="rb")
#   }
#
#   close(f)
#
#   if(verbose) {
#      cat("Best solution at Epoch", epoch.best - 1,
#	 "loss:", losses[epoch.best], "\n")
#   }
#
#   structure(b.best, class="sgd", model=model, lambda1=lambda1,
#	 lambda2=lambda2, features=features, subset=subset,
#	 p=p, epochs=epoch-1, nsamples=length(y), source="file",
#	 losses=losses[-1][2:epoch-1], stepsize=stepsize, anneal=anneal)
#}

rfesgd.character <- function(x, y, p,
      nfeats=sort(unique(c(2^(0:floor(log2(p))), p))), ...)
{
   models <- vector("list", length(nfeats))
   nfeats <- sort(nfeats, decreasing=TRUE)

   for(i in seq(along=nfeats))
   {
      cat(">>>>> nfeats:", nfeats[i], "\n")
      nf <- nfeats[i]

      # If we've already trained a previous model, use it, otherwise start
      # from all features
      n <- ifelse(i == 1, p, nfeats[i-1])

      while(n >= nf)
      {
	 cat("\tTraining", nfeats[i], "\n")
	 # Order features from previous model but ignore intercept
	 r <- if(i > 1) {
	    order(abs(models[[i-1]])[-1], decreasing=TRUE)
	 } else {
	    1:p
	 }
	 cat("\torder:", r, "\n")

	 models[[i]] <- sgd(x, y, p, features=sort(r[1:n]), ...)
	 attr(models[[i]], "order") <- r
	 n <- if(n - nf <= 100) {
	    n - 1
	 } else {
	    n - nf - 100
	 }
      }
   }

   models
}

# Batch gradient descent
gdsvd <- function(x, maxiter=100)
{
   stop("broken")

   # Initialise to anything but zero
   V <- matrix(rnorm(n * m), n, m)
   U <- matrix(rnorm(m * n), m, n)
   stepsize <- 1e-6
   loss <- numeric(N)
   iloss <- numeric(N)
   #loss[1] <- sum((U %*% t(V) - x)^2)
   #iloss[1] <- invloss(U, V, x)
   
   for(i in 1:maxiter)
   {
      #z1 <- 2 * (V %*% t(U) %*% x - diag(n)) %*% t(x) %*% U
      #z2 <- 2 * (V %*% t(U) %*% x - diag(n)) %*% t(V) %*% t(x)
      V <- V - stepsize * z1
      U <- U - stepsize * z2
   
      #xhat <- U %*% t(V) 
      #loss[i] <- sum((xhat - x)^2)
      #iloss[i] <- invloss(U, V, x)
   }

   d <- apply(U, 2, crossprod)
   u <- apply(U, 2, function(x) x / crossprod(x))
   list(u=u, d=d, v=V)
}

# Stochastic gradient descent, looping over samples again when we run out of
# them, total maxiter loops
sgdsvd <- function(x, maxiter=100)
{
   stop("broken")

   m <- nrow(x)
   n <- ncol(x)
   N <- m
   U <- array(rnorm(m * n * N), dim=c(m, n, N))
   V <- array(rnorm(n * n * N), dim=c(n, n, N))
   stepsize <- 1e-3
   loss <- numeric(N)
   loss[1] <- sum((V[, ,1] %*% t(U[, ,1]) - t(x))^2)
   
   for(j in 1:maxiter)
   {
      for(i in 2:N)
      {
         for(k in 1:m)
         {
	    z1 <- U[k, ,i-1] %*% t(V[k, ,i-1])
	    z2 <- V[k, ,i-1] %*% t(U[k, ,i-1])
	    U[k, ,i] <- U[k, ,i-1] - stepsize * 2 * (z1 - x) %*% V[k,,i-1] 
            V[k, ,i] <- V[k, ,i-1] - stepsize * 2 * (z2 - t(x)) %*% U[k,,i-1] 
            xhat <- U[k, ,i-1] %*% t(V[k, ,i-1]) 
            loss[i] <- sum((xhat - x)^2)
            cat("loss:", loss[i], ", cor:", cor(xhat[,1], x[,1]), "\r")
         }
      }
   }


   list(u=u, d=d, v=v)
}

################################################################################
# 
# Utility functions
#

setGeneric("scale", function(object, center, scale, ...) standardGeneric("scale"))

# Get mean and sd from disk file, for each column
# http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
setMethod("scale",
   signature(object="gmatrix", center="ANY", scale="ANY"),
   function(object, center=TRUE, scale=TRUE, verbose=TRUE) {
      mean <- numeric(g@ncol)
      sumsq <- numeric(g@ncol)
      n <- 0
      while(n < g@nrow)
      {
         x <- nextRow(g, loop=FALSE)[[1]]
         m <- length(x) / g@ncol
         if(verbose)
	    cat("Read", n, "row/s\r")
         n <- n + m
         delta <- x - mean
         mean <- mean + delta / n
         sumsq <- sumsq + delta * (x - mean)
      }
      cat("\n")
      if(verbose)
         cat("Read", n, "rows in total\n")
      list(mean=mean, sd=sqrt(sumsq / (n - 1)))
   }
)

test.scale <- function()
{
   set.seed(487134919)
   n <- 100L
   p <- 10L
   x <- matrix(rnorm(n * p), n, p)

   g <- gmatrixMem(x, nrow=n, ncol=p)
   mean <- colMeans(x)
   sd <- apply(x, 2, sd)
   s <- scale(g)

   mean((drop(s$mean) - mean)^2)
   mean((drop(s$sd) - sd)^2)
}

#crossval <- function(nfolds=3, nreps=1, ...)
#{
#
#}

################################################################################
# 
# Prediction functions
#

#setMethod("predict", signature(object="sgd", x="matrix"),
#   function(object, x, scale=TRUE, type=c("linear", "response"), subset=NULL) {
#      
#   }
#)

#setGeneric("predict", function(object, x, ...) standardGeneric("predict"))

setMethod("predict", signature(object="sgd"),
   function(object, g, scale=list(mean=0, norm=1), type=c("linear", "response"), subset=NULL) {

      if(!is(g, "gmatrix"))
	 stop("can only handle x of class gmatrix")

      type <- match.arg(type)
      model <- object@model
      
      if(!model %in% c("logistic", "multinomial") && type == "response")
      {
         stop("don't know what to do with model of type '", attr(b, "model"),
               "and type=response")
      }

      if(is.null(subset))
         subset <- rep(TRUE, g@nrow)

      features <- object@features
      B <- coef(object)
      pr <- matrix(0, g@nrow, ncol(B))

      for(i in 1:g@nrow)
      {
	 cat("reading sample", i, "\r")
         r <- nextRow(g, loop=FALSE)
         if(length(r) == 0)
            stop("Read zero-length data from gmatrix")
         
          if(subset[i]) {
               x <- r[[1]]
               y <- r[[2]]$y
               x <- cbind(1, (x - scale$mean) / scale$norm)[, c(1, features + 1), drop=FALSE]
               pr[i, ] <- x %*% B
          } else if(verbose) {
            cat("skipping", i, "\n")
          }
      }
      cat("\n")
      
      if(type == "linear") {
         pr
      } else if(model == "logistic") {
         1 / (1 + exp(-pr))
      } else if(model == "multinomial") {
         exp(pr) / rowSums(exp(pr))
      }
   }
)

#predict.gd <- function(b, x, scale=TRUE, type=c("linear", "response"),
#      subset=NULL)
#{
#   type <- match.arg(type)
#   if(!attr(b, "model") %in% c("logistic", "multinomial")
#	 && type == "response")
#   {
#      stop("don't know what to do with model of type '",
#	 attr(b, "model"), "' and type='response'")
#   }
#
#   if(is.null(subset))
#      subset <- rep(TRUE, nrow(x))
#
#   if(scale)
#      x <- scale(x)
#   pr <- cbind(1, x[subset, , drop=FALSE]) %*% b
#   if(type == "linear") {
#      pr
#   } else if(attr(b, "model") == "logistic") {
#      1 / (1 + exp(-pr))
#   } else if(attr(b, "model") == "multinomial") {
#      exp(pr) / rowSums(exp(pr))
#   }
#}
#
#predict.sgd <- function(b, x, blocksize=1, type=c("linear", "response"),
#      subset=NULL, verbose=TRUE)
#{
#   type <- match.arg(type)
#   if(is.matrix(x))
#      return(predict.gd(b, x, scale=FALSE, type=type, subset=subset))
#
#   if(nchar(x) == 0)
#      stop("filename x not supplied")
#   
#   if(!attr(b, "model") %in% c("logistic", "multinomial")
#      && type == "response")
#   {
#      stop("don't know what to do with model of type '", attr(b, "model"),
#	    "and type=response")
#   }
#
#   fname <- x
#   f <- file(fname, open="rb")
#   pr <- numeric(0)
#   p <- attr(b, "p")
#   if(is.null(subset))
#      subset <- rep(TRUE, attr(b, "nsamples"))
#
#   features <- attr(b, "features")
#
#   i <- 1
#   while(TRUE)
#   {
#      dat <- readBin(f, what="numeric", n=p)
#      if(length(dat) == 0)
#	 break
#      
#      if(subset[i]) {
#         x <- matrix(dat, nrow=blocksize, ncol=p)[, features, drop=FALSE]
#         pr <- c(pr, drop(cbind(1, x) %*% b))
#      } else if(verbose) {
#         cat("skipping", i, "\n")
#      }
#      
#      i <- i + 1
#   }
#   close(f)
#   
#   if(type == "linear") {
#      pr
#   } else if(attr(b, "model") == "logistic") {
#      1 / (1 + exp(-pr))
#   } else if(attr(b, "model") == "multinomial") {
#      exp(pr) / rowSums(exp(pr))
#   }
#}
#



#gd <- function(x, y, loss, dloss, threshold=1e-4, stepsize=1 / nrow(x),
#      maxiter=1000, initial=rep(1, ncol(x)))
#{
#   p <- ncol(x)
#   b <- matrix(0, maxiter, p)
#   b[1, ] <- initial
#   l.old <- loss(x, y, b[1, ])
#   l.new <- Inf
#   i <- 1
#
#   while(i < maxiter && any(is.finite(c(l.old, l.new)))
#	 && max(abs(l.old - l.new)) >= threshold)
#   {
#      l.old <- l.new
#      b[i + 1, ] <- drop(b[i, ]  - stepsize * dloss(x, y, b[i, ]))
#      cat(b[i + 1, ], "\r")
#      l.new <- loss(x, y, b[i + 1, ])
#      i <- i + 1
#   }
#   b[1:i, ]
#}

#sgd <- function(x, y, loss, dloss, threshold=1e-4, stepsize=1 / nrow(x),
#      initial=rep(1, ncol(x)))
#{
#   p <- ncol(x)
#   n <- nrow(x)
#   b <- matrix(0, n, p)
#   b[1, ] <- initial
#
#   #l.old <- l2loss(x, y, b)
#   #l.new <- Inf
#   #while(max(abs(l.old - l.new)) >= threshold)
#   for(i in 1:1000)
#   {
#      #l.old <- l.new
#      for(j in 1:(n-1))
#      {
#	 b[j + 1, ] <- drop(b[j, ] - stepsize * dloss(x[j, , drop=FALSE], y[j], b[j, ]))
#      }
#      #l.new <- l2loss(x, y, b)
#      #i <- i + 1
#   }
#   b
#}


