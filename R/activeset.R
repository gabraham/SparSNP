library(glmnet)

options(error=dump.frames)

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

# assumes scaled inputs
step.linear <- function(xj, lp, y)
{
   drop(crossprod(xj, lp - y)) / (length(y) - 1)
}

# assumes scaled inputs
step.logis <- function(xj, lp, y)
{
   p <- 1 / (1 + exp(-lp))
   d1 <- mean(xj * (p - y))
   q <- p * (1 - p)
   d2 <- mean(xj^2 * q)
   d1 / d2
}

step.logis2 <- function(xj, lp, y)
{
   
}

# assumes scaled inputs
step.sqrhinge <- function(xj, lp, y)
{
   z <- y * lp - 1

   mean(y * xj * z * (z < 0))
}

loss.linear <- function(lp, y)
{
   lp <- cbind(lp)
   apply(lp, 2, function(l) {
      mean((l - y)^2)
   })
}

loss.logis <- function(lp, y)
{
   lp <- cbind(lp)
   apply(lp, 2, function(l) {
      mean(log(1 + exp(l * y)) - y * l)
   })
}

loss.sqrhinge <- function(lp, y)
{
   mean(pmax(0, 1 - y * lp)^2)
}

cd.plain <- function(x, y, step=step.linear, loss=loss.linear,
   lambda1=0, lambda1max=0,
   beta=numeric(ncol(x)), lp=numeric(nrow(x)),
   maxepoch=50, maxiter=50, eps=1e-4)
{
   x <- cbind(1, x)
   p <- ncol(x)
   n <- nrow(x)
   loss_new <- loss_old <- Inf

   grad <- matrix(0, maxepoch, p)
   LP <- matrix(0, maxepoch, n)
   
   for(epoch in 1:maxepoch)
   {
      for(j in 1:p)
      {
	 for(iter in 1:maxiter)
	 {
	    s <- step(x[,j], lp, y)
	    beta_new <- beta[j] - s
	    if(j > 1)
	       beta_new <- soft(beta_new, lambda1)

	    diff <- beta_new - beta[j]
	    lp <- lp + x[, j] * diff
	    beta[j] <- beta_new
	    loss_new <- loss(lp, y)
	    
	    if(abs(loss_new - loss_old) < 1e-2 || abs(diff) < eps)
	       break
	    loss_old <- loss_new
	 }
      }
      grad[epoch, ] <- crossprod(x, lp-y)
      LP[epoch, ] <- lp
   }
   list(beta=beta, grad=grad, LP=LP)
}

cd.activeset <- function(x, y, step, loss, lambda1=0, beta=numeric(ncol(x)),
   lp=numeric(nrow(x)), maxepoch=1000, maxiter=100, eps=1e-7, intercept=TRUE)
{
   p <- ncol(x)
   n <- nrow(x)
   beta_old <- beta
   active_new <-  rep(TRUE, p)
   active_old <- rep(FALSE, p)
   allconverged <- 0
   loss_old <- 1e9
   loss_new <- 1e10
   loss_null <- loss(numeric(n) + mean(y), y)

   cat("lambda1:", lambda1, "\n")
   for(epoch in 1:maxepoch)
   {
      for(j in 1:p)
      {
	 if(active_new[j]) 
	 {
	    conv <- FALSE
	    for(iter in 1:maxiter)
	    {
	       s <- step(x[,j], lp, y)
               beta_new <- beta[j] - s
               if(intercept && j > 1)
                  beta_new <- soft(beta_new, lambda1)

	       diff <- beta_new - beta[j]
	       #lp <- pmin(10, pmax(lp + x[, j] * diff, -10))
	       lp <- lp + x[,j] * diff
	       beta[j] <- beta_new
	       #loss_new <- loss(lp, y)
	       #loss_ratio <- loss_new / loss_null
	       #if(abs(diff) < eps 
	       #	  || abs(loss_old - loss_new) <= 1e-3
	       #	  || loss_new <= 1e-3
	       #	  || loss_ratio >= 0.999)
	          break
	    }
	    if(iter == maxiter)
	       cat("maxiter (", maxiter, ") reached for variable", j, "\n")
	    active_new[j] <- beta[j] != 0
	 }
      }

      if(all(abs(beta - beta_old) <= eps))
      {
         allconverged <- allconverged + 1
	 cat("allconverged:\n", allconverged)

         if(allconverged == 1) {
	    cat("allconverged 1\n")
            active_old <- active_new
            active_new[] <- TRUE
         } else if(allconverged == 2) {
	    cat("allconverged 2\n")
            if(all(active_old == active_new))
	    {
	       cat("converged\n")
               break
	    }
            active_old <- active_new
            active_new[] <- TRUE
            allconverged <- 1
         }
      }
      else
      {
	 cat(epoch, "allconverged 0\n")
	 allconverged <- 0 
      }

      if(epoch == maxepoch)
      {
	 cat("failed to converged within", maxepoch, "epochs\n")
	 browser()
      }

      beta_old <- beta
      loss_old <- loss_new
   }

  
   cat("\n")
   beta
}

cd.activeset.experimental <- function(x, y, step, loss, lambda1=0,
   beta=numeric(ncol(x)), lp=numeric(nrow(x)), maxepoch=500,
   maxiter=50, eps=1e-10)
{
   p <- ncol(x)
   n <- nrow(x)
   active_new <-  rep(TRUE, p)
   active_old <- rep(FALSE, p)
   allconverged <- 0
   numconverged <- 0

   cat("lambda1:", lambda1, "\n")
   for(epoch in 1:maxepoch)
   {
      numconverged = 0
      for(j in 1:p)
      {
	 conv <- TRUE
	 if(active_new[j]) 
	 {
	    conv <- FALSE
	    for(iter in 1:maxiter)
	    {
	       s <- step(x[,j], lp, y)
               beta_new <- beta[j] - s
               if(j > 1)
                  beta_new <- soft(beta_new, lambda1)

	       s <- beta_new - beta[j]
	       lp <- lp + x[, j] * s
	       beta[j] <- beta_new

	       if(abs(s) < eps)
	       {
		  conv <- TRUE
	          break
	       }
	    }
	    if(iter == maxiter)
	       cat("maxiter reached for variable", j, "\n")
	    active_new[j] <- beta[j] != 0
	 }
	 numconverged <- numconverged + conv
      }

      if(numconverged == p)
      {
         allconverged <- allconverged + 1

         if(allconverged == 1) {
            active_old <- active_new
            #active_new[] <- TRUE
         } else if(allconverged == 2) {
            if(all(active_old == active_new))
               break
            active_old <- active_new
            #active_new[] <- TRUE
            allconverged <- 1
         }
      }
      else
	 allconverged <- 0 
   }

   if(epoch == maxepoch)
      cat("failed to converged within", maxepoch, "epochs\n")

   cat("\n")
   beta
}

cd.activeset.broken <- function(x, y, step, loss, lambda1=0, beta=numeric(ncol(x)),
   lp=numeric(nrow(x)), maxepoch=500, maxiter=50, eps=1e-4)
{
   p <- ncol(x)
   n <- nrow(x)
   beta_old <- beta
   active <-  rep(TRUE, p)
   conv <-  rep(TRUE, p)
   allconverged <- 0
   numconverged <- 0

   cat("lambda1:", lambda1, "\n")
   for(epoch in 1:maxepoch)
   {
      numconverged <- 0
      for(j in 1:p)
      {
	 if(active[j]) 
	 {
	    for(iter in 1:maxiter)
	    {
	       s <- step(x[,j], lp, y)
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
	    active[j] <- beta[j] != 0
	 }
	 conv[j] <- abs(beta[j] - beta_old[j]) <= eps
	 beta_old[j] <- beta[j]
	 numconverged <- numconverged + conv[j]
      }

      if(numconverged == p)
      {
	 allconverged <- allconverged + 1
	 if(allconverged == 1)
	 {
	    active[] <= TRUE
	 }
	 else
	 {
	    break;
	 }
      }
      else
      {
	 allconverged <- 0
	 active[] <- TRUE
      }
   }

   if(epoch == maxepoch)
      cat("failed to converged within", maxepoch, "epochs\n")
   cat("\n")

   cat("finished in", epoch, "epochs with", sum(beta != 0), "active vars\n")
   beta
}

# Slightly different formulation
cd.activeset2 <- function(x, y, step, loss, lambda1=0, lambda1max=0,
   beta=numeric(ncol(x)), lp=numeric(nrow(x)),
   maxepoch=1000, maxiter=100, eps=1e-7)
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
	       #s <- sum(x[,j] * (lp - y)) / n
	       s <- step(x[,j], lp, y)
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


cd.activeset.20110921 <- function(x, y, step, loss, lambda1=0,
   beta=numeric(ncol(x)), lp=numeric(nrow(x)),
   maxepoch=1000, maxiter=50, eps=1e-4, intercept=TRUE)
{
   p <- ncol(x)
   n <- nrow(x)
   beta_old <- beta
   active_new <-  rep(TRUE, p)
   active_old <- rep(FALSE, p)
   allconverged <- 0
   loss_old <- 1e9
   loss_new <- 1e10
   loss_null <- loss(numeric(n) + mean(y), y)

   lp_old <- lp
   # check unpenalised marginal gradient
   grad_null <- numeric(p)
   #for(j in 1:p)
   #   grad_null[j] <- abs(step(x[,j], lp_old, y)) > lambda1
   S <- B <- matrix(0, p, maxepoch)

   cat("lambda1:", lambda1, "\n")
   for(epoch in 1:maxepoch)
   {
      discord <- 0
      for(j in 1:p)
      {
	 if(active_new[j])
	 {
	    conv <- FALSE
	    for(iter in 1:maxiter)
	    {
	       s <- step(x[,j], lp, y)
	       S[j, epoch] <- s
               beta_new <- beta[j] - s
               if(intercept && j > 1)
                  beta_new <- soft(beta_new, lambda1)

	       diff <- beta_new - beta[j]
	       #lp <- pmin(10, pmax(lp + x[, j] * diff, -10))
	       lp <- lp + x[, j] * diff
	       beta[j] <- beta_new
	       B[j, epoch] <- beta_new
	#       loss_new <- loss(lp, y)
	#       loss_ratio <- loss_new / loss_null
	#       if(any(is.nan(c(loss_new, loss_new))))
	#	  browser()
	#       if(abs(diff) < eps 
	#	  || abs(loss_old - loss_new) <= 1e-2
	#	  || loss_new <= 1e-2
	#	  || loss_ratio >= 0.99)
	          break
	    }
	    #if(iter == maxiter)
	     #  cat("maxiter (", maxiter, ") reached for variable", j, "\n")
	    active_new[j] <- beta[j] != 0
	    #d <- as.integer(active_new[j] != grad_null[j]);
	    #if(d != 0)
	    #   browser()
	    #discord <- discord + d
	 }
      }

      if(all(abs(beta - beta_old) <= eps))
      {
         allconverged <- allconverged + 1
	 cat("allconverged:\n", allconverged)

         if(allconverged == 1) {
	    cat("allconverged 1\n")
            active_old <- active_new
            active_new[] <- TRUE
         } else if(allconverged == 2) {
	    cat("allconverged 2\n")
            if(all(active_old == active_new))
	    {
	       cat("converged\n")
               break
	    }
            active_old <- active_new
            active_new[] <- TRUE
            allconverged <- 1
         }
      }
      else
      {
	 cat(epoch, "allconverged 0\n")
	 allconverged <- 0 
      }

      if(epoch == maxepoch)
      {
	 cat("failed to converged within", maxepoch, "epochs\n")
	 browser()
      }

      beta_old <- beta
      loss_old <- loss_new
   }

   cat("epoch:", epoch, "\n")
   S <- S[, 1:epoch]
   B <- B[, 1:epoch]
browser()
  
   cat("\n")
   #cat("discord:", discord, "\n")
   beta
}

cd.activeset.20110922 <- function(x, y, step, loss, lambda1=0,
   beta=numeric(ncol(x)), lp=numeric(nrow(x)),
   maxepoch=1000, maxiter=100, eps=1e-7, intercept=TRUE)
{
   p <- ncol(x)
   n <- nrow(x)
   beta_old <- beta
   active_new <-  rep(TRUE, p)
   active_old <- rep(FALSE, p)
   allconverged <- 0
   loss_old <- 1e9
   loss_new <- 1e10
   loss_null <- loss(numeric(n) + mean(y), y)

   # check unpenalised marginal gradient
   grad_null <- numeric(p)
   for(j in 1:p)
      grad_null[j] <- abs(step(x[,j], lp, y)) > lambda1
   if(intercept)
      grad_null[j] <- TRUE

   cat("lambda1:", lambda1, "\n")
   for(epoch in 1:maxepoch)
   {
      for(j in 1:p)
      {
	 if(active_new[j] && grad_null[j])
	 {
	    conv <- FALSE
	    for(iter in 1:maxiter)
	    {
	       s <- step(x[,j], lp, y)
               beta_new <- beta[j] - s
               if(intercept && j > 1)
                  beta_new <- soft(beta_new, lambda1)

	       diff <- beta_new - beta[j]
	       #lp <- pmin(10, pmax(lp + x[, j] * diff, -10))
	       lp <- lp + x[,j] * diff
	       beta[j] <- beta_new
	       #loss_new <- loss(lp, y)
	       #loss_ratio <- loss_new / loss_null
	       #if(any(is.nan(c(loss_new, loss_new))))
	       #   browser()
	       #if(abs(diff) < eps 
	       #   || abs(loss_old - loss_new) <= 1e-3
	       #   || loss_new <= 1e-3
	       #   || loss_ratio >= 0.999)
	       #   break
	       break
	    }
	    if(iter == maxiter)
	       cat("maxiter (", maxiter, ") reached for variable", j, "\n")
	    active_new[j] <- beta[j] != 0
	 }
      }

      if(all(abs(beta - beta_old) <= eps))
      {
         allconverged <- allconverged + 1
	 cat("allconverged:\n", allconverged)

         if(allconverged == 1) {
	    cat("allconverged 1\n")
            active_old <- active_new
            active_new[] <- TRUE
         } else if(allconverged == 2) {
	    cat("allconverged 2\n")
            if(all(active_old == active_new))
	    {
	       cat("converged\n")
               break
	    }
            active_old <- active_new
            active_new[] <- TRUE
            allconverged <- 1
         }
      }
      else
      {
	 cat(epoch, "allconverged 0\n")
	 allconverged <- 0 
      }

      if(epoch == maxepoch)
      {
	 cat("failed to converged within", maxepoch, "epochs\n")
	 browser()
      }

      beta_old <- beta
      loss_old <- loss_new
   }

  
   cat("\n")
   beta
}

cd.group.plain <- function(x, y, step, loss, lambda1=0, lambda1max=0,
   beta=numeric(ncol(x)), lp=numeric(nrow(x)),
   maxepoch=500, maxiter=50, eps=1e-4)
{
   p <- ncol(x)
   n <- nrow(x)
   #lp <- numeric(n)
   #beta <- numeric(p)
   loss_new <- loss_old <- Inf
   
   for(epoch in 1:maxepoch)
   {
      for(j in 1:p)
      {
	 for(iter in 1:maxiter)
	 {
	    s <- step(x[,j], lp, y)
	    beta_new <- beta[j] - s
	    if(j > 1)
	       beta_new <- soft(beta_new, lambda1)

	    diff <- beta_new - beta[j]
	    lp <- lp + x[, j] * diff
	    beta[j] <- beta_new
	    loss_new <- loss(lp, y)
	    if(abs(loss_new - loss_old) < 1e-2 || abs(diff) < eps)
	       break
	    loss_old <- loss_new
	 }
      }
   }
   beta
}


#check.kkt <- function(x, y, beta, step, lambda)
#{
#   s <- numeric(ncol(x))
#   for(j in 1:p)
#   {
#      s[j] <- abs(step(x, y)) 
#   }
#}

