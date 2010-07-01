# Test CD

library(glmnet)

set.seed(43249210)

n <- 1000
p <- 10000
w <- sample(0:1, p + 1, replace=TRUE, prob=c(0.9, 0.1)) 
beta <- rnorm(p + 1) * w

x <- matrix(sample(c(0, 1, 2), size=n * p, replace=TRUE), n, p)
y <- ifelse(runif(n) <= plogis(cbind(1, x) %*% beta + rnorm(n, 1, 2)), 1, 0)
table(y)

z <- cbind(y, x)
writeBin(as.raw(z), con="x.bin.t")


base::q("no")

source("../tests/optim.R")


#b1 <- cd1(x, y, logloss, logdloss, logd2loss)
#b2 <- cd2(x, y, logloss, logd1phi, logd2phi)

d1phi <- logd1phi
d2phi <- logd2phi
g2 <- glmnet(x, factor(y), family="binomial") 

getlambda1max <- function(x, y)
{
   x <- cbind(1, x)
   n <- nrow(x)
   p <- ncol(x)
   lp <- x %*% beta
   
   s <- numeric(p)
   for(j in 1:p)
   {
      grad <- sum(x[, j] * (d1phi(lp) - y))
      d2 <- sum(x[,j]^2 * d2phi(lp))
      
      beta.old <- beta[j]
      d <- beta.old - grad / d2
      cat(j, grad, d2, d, "\n")
      if(j == 1) {
	 beta[j] <- d
      }
   
      lp <- lp + x[, j] * (beta[j] - beta.old)
      s[j] <- d
   }
   
   max(abs(s[-1]))
}

l1max <- getlambda1max(x, y)

nl1 <- 50
l1min <- 1e-3 * l1max
s <- (log(l1max) - log(l1min)) / nl1
l <- exp(log(l1max) - s * 0:nl1)
b21 <- sapply(l, cd2, x=xs, y=y, lossfunc=logloss, d1phi=logd1phi, d2phi=logd2phi)

df <- apply(b21[-1,], 2, function(x) sum(x != 0))
#l1max2 <- l[which(c(0, diff(df)) > 0)]

plot(l, df, log="x")

stop()

nl1 <- 50
cmd <- sprintf("time ./cd -model logistic -f x.bin.t -n %s -p %s -beta beta.csv \\
-epochs 100 -thresh 1e-4 -v -nl1 %s", n, p, nl1)
system(cmd)

b.cd <- sapply(1:nl1 - 1, function(i) {
   read.csv(sprintf("beta.csv.%s", i), header=FALSE)[,1]
})
b.cd.z <- apply(b.cd[-1,], 2, function(x) sum(x == 0))
b.cd.1 <- b.cd[, which.min((b.cd.z - sum(w == 0))^2)]

g1 <- glm(y ~ x, family=binomial)
b.glm <- coef(g1)
b.glmnet <- as.matrix(coef(g2))
b.glmnet.z <- apply(b.glmnet[-1,], 2, function(x) sum(x == 0))
b.glmnet.1 <- b.glmnet[, which.min((b.glmnet.z - sum(w == 0))^2)]

#b.cd2 <- cd2(x, y, logloss, logd1phi, logd2phi, maxiter=100)

res <- cbind(beta=beta, glm=b.glm, cd=b.cd.1, glmnet=b.glmnet.1)
pairs(res, upper.panel=function(x, y, ...) {
   points(x, y, ...)
   abline(0, 1)
})
cor(res)

matplot(res, type="l", lty=1)

(rmse.beta <- apply(res, 2, function(x) sqrt(mean(res[,1] - x)^2)))
(rmse.glm <- apply(res, 2, function(x) sqrt(mean(res[,2] - x)^2)))

