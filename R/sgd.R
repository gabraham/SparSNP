# Squared loss for regression
l2loss <- function(x, y, b)
{
   sum((y - x %*% b)^2)
}

l2dloss <- function(x, y, b)
{
   crossprod(y - x %*% b, -x)
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

logdloss <- function(x, y, b)
{
   p <- exp(x %*% b) 
   crossprod(-x, y - p / (1 + p))
}

# Multinomial logistic
# y is an indicator matrix
# b is a matrix of coefficients
mlogloss <- function(x, y, b)
{
   
}

mlogdloss <- function(x, y, b)
{
   p <- exp(x %*% b)
   crossprod(-x, y - p / rowSums(p))
}

# Cox proportional hazards
coxloss <- function(x, y, b)
{
   
}

coxdloss <- function(x, y, b)
{
}

# Batch gradient descent
gd <- function(x, y, model=c("linear", "logistic", "hinge"),
      lambda1=0, lambda2=0, lambdaE=0, alpha=NULL,
      threshold=1e-4, stepsize=1 / nrow(x), anneal=1e-9, maxiter=100,
      scale=TRUE)
{
   model <- match.arg(model)
   loss <- switch(model, linear=l2loss, logistic=logloss, hinge=hingeloss)
   dloss <- switch(model, linear=l2dloss, logistic=logdloss, hinge=hingedloss)

   x <- if(scale) {
      cbind(1, scale(x))
   } else cbind(1, x)
   p <- ncol(x)
   b <- rep(0, p)
   l.old <- Inf
   l.new <- 0

   # Elastic net
   if(lambdaE > 0 && !is.null(alpha))
   {
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
      grad <- drop(dloss(x, y, b)) + lambda2 * b + lambda1 * sign(b)
      b <- b - stepsize * grad
      l.new <- loss(x, y, b)
      cat(l.new, "\n")
      stepsize <- stepsize * 1 / (1 + anneal)
      i <- i + 1
   }
   cat("converged at", i-1, "iterations\n")

   structure(b, class="gd", model=model, iter=i)
}

sgd <- function(x, y, ...)
{
   UseMethod("sgd")
}

sgd.matrix <- function(x, y, model=c("linear", "logistic", "hinge"),
      lambda1=0, lambda2=0, lambdaE=0,
      alpha=NULL, threshold=1e-4, stepsize=1 / nrow(x),
      maxepochs=1, anneal=1e-9, maxiter=Inf, scale=TRUE)
{
   model <- match.arg(model)
   loss <- switch(model, linear=l2loss, logistic=logloss, hinge=hingeloss)
   dloss <- switch(model, linear=l2dloss, logistic=logdloss, hinge=hingedloss)

   x <- if(scale) {
      cbind(1, scale(x))
   } else {
      cbind(1, x)
   }
   p <- ncol(x)
   n <- nrow(x)
   n <- min(maxiter, n)
   b <- rep(0, p)

   # Don't penalise intercept
   lambda1 <- c(0, rep(lambda1, p - 1))
   lambda2 <- c(0, rep(lambda2, p - 1))

   l.old <- Inf
   l.new <- 0
   epoch <- 1
   while(epoch <= maxepochs && max(abs(l.old - l.new)) >= threshold)
   {
      l.new <- l.new <- 0
      for(i in 1:n)
      {
	 grad <- (drop(dloss(x[i, , drop=FALSE], y[i], b)) 
	    + lambda2 * b + lambda1 * sign(b))
	 l.new <- l.new + loss(x[i, , drop=FALSE], y[i], b)
	 b <- b - stepsize * grad
      }
      stepsize <- stepsize * 1 / (1 + anneal)
      cat("Epoch", epoch, ", loss:", l.new, "\n")
      epoch <- epoch + 1
   }
   structure(b, class="sgd", model=model, p=p, epochs=epoch, source="matrix") 
}

sgd.character <- function(x="", y, p, model=c("linear", "logistic", "hinge"),
      lambda1=0, lambda2=0, lambdaE=0, alpha=NULL, threshold=1e-4,
      stepsize=1e-4, maxepochs=1, anneal=stepsize, blocksize=1,
      maxiter=Inf, subset=NULL, saveloss=FALSE)
{
   if(nchar(x) == 0)
      stop("filename x not supplied")

   if(!is.null(subset) && blocksize > 1)
      stop("blocksize > 1 currently not supported when subset is used")

   if(is.null(subset))
      subset <- rep(TRUE, length(y))

   model <- match.arg(model)
   loss <- switch(model, linear=l2loss, logistic=logloss, hinge=hingeloss)
   dloss <- switch(model, linear=l2dloss, logistic=logdloss, hinge=hingedloss)

   fname <- x
   f <- file(fname, open="rb")

   # Don't penalise intercept
   lambda1 <- c(0, rep(lambda1, p))
   lambda2 <- c(0, rep(lambda2, p))

   b <- rep(0, p + 1)
   losses <- rep(0, maxepochs + 1)
   totalloss <- 0
   epoch <- 2
   while(epoch <= 2 || (epoch <= maxepochs + 1
	 && max(abs(losses[epoch] - losses[epoch-1])) >= threshold))
   {
      #losses[epoch] <- losses[epoch-1]
      i <- 1
      while(TRUE)
      {
	 #cat("block", i, "\r")
	 dat <- readBin(f, what="numeric", n=p * blocksize)
	 if(length(dat) == 0 || i > maxiter)
	    break

	 if(subset[i]) {
	    x <- matrix(dat, nrow=min(blocksize, length(dat) / p), ncol=p,
	          byrow=TRUE)
	    x <- cbind(1, x)
	    k <- ((i-1) * blocksize + 1):(i * blocksize)
	    k <- k[1:nrow(x)]
	    l <- loss(x, y[k], b)
	    losses[epoch] <- losses[epoch] + l
	    if(l > 0)
	       cat(i, l, "\n")
	    grad <- drop(dloss(x, y[k], b)) + lambda2 * b + lambda1 * sign(b)
	    b <- b - stepsize * grad
	 } else {
	    cat("skipping", i, "\n")
	 }
	 i <- i + 1
      }
      stepsize <- stepsize / (1 + anneal)
      cat("Epoch", epoch-1, ", loss:", losses[epoch], "\n")
      epoch <- epoch + 1
      close(f)
      f <- file(fname, open="rb")
   }

   close(f)

   structure(b, class="sgd", model=model, p=p, epochs=epoch-1, source="file",
	 losses=losses[-1], stepsize=stepsize, anneal=anneal)
}

#crossval <- function(nfolds=3, nreps=1, ...)
#{
#
#}

predict.gd <- function(b, x, scale=TRUE, type=c("linear", "response"),
      subset=NULL)
{
   type <- match.arg(type)
   if(attr(b, "model") != "logistic" && type == "response")
   {
      stop("don't know what to do with model of type '", attr(b, "model"),
	    "and type=response")
   }

   if(is.null(subset))
      subset <- rep(TRUE, nrow(x))

   if(scale)
      x <- scale(x)
   pr <- cbind(1, x[subset, , drop=FALSE]) %*% b
   if(type == "linear") {
      pr
   } else if(attr(b, "model") == "logistic") {
      1 / (1 + exp(-pr))
   }
}

predict.sgd <- function(b, x, blocksize=1, type=c("linear", "response"),
      subset=NULL)
{
   type <- match.arg(type)
   if(is.matrix(x))
      return(predict.gd(b, x, scale=FALSE, type=type, subset=subset))

   if(nchar(x) == 0)
      stop("filename x not supplied")
   
   if(attr(b, "model") != "logistic" && type == "response")
   {
      stop("don't know what to do with model of type '", attr(b, "model"),
	    "and type=response")
   }

   fname <- x
   f <- file(fname, open="rb")
   pr <- numeric(0)
   p <- attr(b, "p")
   if(is.null(subset))
      subset <- rep(TRUE, p)

   i <- 1
   while(TRUE)
   {
      dat <- readBin(f, what="numeric", n=p)
      if(length(dat) == 0)
	 break
      
      if(subset[i]) {
         x <- matrix(dat, nrow=blocksize, ncol=p) 
         pr <- c(pr, drop(cbind(1, x) %*% b))
      } else {
         cat("skipping", i, "\n")
      }
      
      i <- i + 1
   }
   close(f)
   
   if(type == "linear") {
      pr
   } else if(attr(b, "model") == "logistic") {
      1 / (1 + exp(-pr))
   }
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


