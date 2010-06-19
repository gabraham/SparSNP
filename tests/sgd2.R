# SGD proving ground

library(ElemStatLearn)
library(glmnet)

logloss <- function(x, y, b)
{
   p <- x %*% b
   log(1 + exp(p)) - y * p
}

logdloss <- function(x, y, b)
{
   p <- exp(x %*% b)
   x %*% (p / (1 + p) - y)
}

l2loss <- function(x, y, b)
{
   sum((y - x %*% b)^2)
}

l2dloss <- function(x, y, b)
{
   x * (x %*% b - y)
}

threshold <- function(a, b)
{
   x <- a - b
   x[sign(a) != 0 & sign(a) != sign(x)] <- 0
   x
}

truncate <- function(x, epsilon=0)
{
   x[abs(x) <= epsilon] <- 0
   x
}

sgd <- function(x, y, lossfunc=l2loss, dlossfunc=l2dloss,
   updatefunc=1, maxepochs=100, stepsize=1e-5, lambda1=0)
{
   x <- cbind(1, x)
   n <- length(y)
   p <- ncol(x)
   b <- numeric(p)
   loss <- loss2 <- rep(as.numeric(NA), maxepochs)
   anneal <- stepsize

   # Plain updates
   update1 <- function(epoch)
   {
      curloss <- 0
      for(i in 1:n)
      {
         curloss <- curloss + lossfunc(x[i,], y[i], b)
         grad <- dlossfunc(x[i,], y[i], b)
	 b <- b - stepsize * drop(grad + lambda1 * sign(b))
      }
      cat(epoch, curloss, "\n")
     
      list(curloss, b, stepsize)
   }

   # Step halving
   update2 <- function(epoch)
   {
      curloss <- 0
      for(i in 1:n)
      {
	 curloss <- curloss + lossfunc(x[i,], y[i], b)
         grad <- dlossfunc(x[i,], y[i], b)
	 b <- b - stepsize * drop(grad + lambda1 * sign(b))
      }
      cat(epoch, curloss, "\n")
      if(epoch > 1 && curloss > loss[epoch-1])
         stepsize <- stepsize / 2 

      list(curloss, b, stepsize)
   }

   # Annealing
   update3 <- function(epoch)
   {
      curloss <- 0
      for(i in 1:n)
      {
	 curloss <- curloss + lossfunc(x[i,], y[i], b)
         grad <- dlossfunc(x[i,], y[i], b)
	 b <- b - stepsize * drop(grad + lambda1 * sign(b))
      }
      cat(epoch, curloss, "\n")
      stepsize <- stepsize / (1 + anneal)

      list(curloss, b, stepsize)
   }

   updatefunc <- switch(updatefunc, update1, update2, update3)
   
   for(epoch in 1:maxepochs)
   {
      u <- updatefunc(epoch)
      loss[epoch] <- u[[1]]
      b <- u[[2]]
      stepsize <- u[[3]]
      loss2[epoch] <- lossfunc(x, y, b)
   }

   list(loss=loss, loss2=loss2, beta=b)
}

y <- prostate$lpsa[prostate$train]
x <- scale(as.matrix(prostate[, 1:8]))[prostate$train, ]
#x <- as.matrix(prostate[, 1:8])[prostate$train, ]
n <- nrow(x)
p <- ncol(x)

s <- sgd(x, y, updatefunc=1, stepsize=1e-5, maxepochs=100)

stop()

g <- glmnet(x, y, family="gaussian")

l <- lm(y ~ x)

#lambda <- c(2^(9:-9), 0)
lambda <- g$lambda
nl <- length(lambda)

updates <- c("plain"=1, "step halving"=2, "annealing"=3)
res <- lapply(updates, function(up) {
   lapply(lambda, function(l) {
      sgd(x, y, lambda1=l, maxepochs=500, stepsize=0.25 / n, updatefunc=up)
   })
})


# SMIDAS
f <- "smidas_x.dat"
cat(dim(x), "\n", file=f)
for(i in 1:n)
{
   cat(p, paste(1:n - 1, x[i,]), "\n", file=f, append=TRUE)
}

write.table(y, file="smidas_y.dat", row.names=FALSE, col.names=FALSE)

res.smd <- lapply(seq(along=lambda), function(i) {
   lapply(1:100, function(j) {
      system(sprintf("~/Software/SMIDAS/smidas -iters %s \\
-printout 1 -loss 1 -lambda %s smidas_x.dat smidas_y.dat > smidas.%s",
	 j, lambda[i], i))
      loss <- read.table(pipe(
	 sprintf("head -%s smidas.%s | tail -1 | cut -f 5 -d ' ' ",
	    j, i)))
      beta <- read.table(pipe(
	 sprintf("tail -%s smidas.%s", p, i)))
      list(loss=loss, beta=c(0, beta))
   })
})

res <- c(res, res.smd)

beta <- lapply(res, function(r) sapply(r, function(x) x$beta)[-1,])
loss <- lapply(seq(along=lambda), function(i) {
   sapply(res, function(r) r[[i]]$loss)
})

ylim <- range(coef(l)[-1], unlist(beta))

pdf("prostate2.pdf", width=12)
par(mfrow=c(1, 3))

for(i in seq(along=beta))
{
   matplot(seq(along=lambda), t(beta[[i]]), type="b", pch=21, cex=0.5,
      xlim=c(1, length(lambda) + 1), ylim=ylim, xlab="Penalty", xaxt="n",
      main=names(beta)[i], ylab="Weight", col=1:p)
   axis(side=1, at=1:length(lambda), labels=lambda)
   points(rep(length(lambda) + 1, p), coef(l)[-1], pch=19, cex=1, col=1:p)
}

par(mfrow=c(1, 3))

for(i in c(1, 10, 20))
{
   matplot(loss[[i]], type="b", pch=21, cex=0.5,
	 lty=1, xlab="Epochs", ylab="Loss",
	 main=paste("lambda=", lambda[i], sep=""))
   par(usr=c(0, 1, 0, 1))
   legend(x=0.65, y=0.99, legend=names(updates), col=1:3, lwd=3)
}

dev.off()

