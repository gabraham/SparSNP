# Test CD

library(glmnet)

set.seed(43249210)

n <- 20
p <- 5
w <- sample(0:1, p + 1, replace=TRUE, prob=c(0.9, 0.1)) 
beta <- rnorm(p + 1) * w

x <- matrix(sample(c(0, 1, 2), size=n * p, replace=TRUE), n, p)
y <- ifelse(runif(n) <= plogis(cbind(1, x) %*% beta + rnorm(n, 1, 2)), 1, 0)
table(y)

z <- cbind(y, x)
writeBin(as.raw(z), con="x.bin.t")


#base::q("no")

source("../tests/optim.R")


#x2 <- scale(x)
x2 <- x
#x2[, apply(x, 2, var) == 0] <- 0

#b1 <- cd1(x, y, logloss, logdloss, logd2loss)
#b2 <- cd2(x, y, logloss, logd1phi, logd2phi)

d1phi <- logd1phi
d2phi <- logd2phi
g2 <- glmnet(x2, factor(y), family="binomial") 
g2$lambda[1]

getlambda1max <- function(x, y)
{
   x <- cbind(1, x)
   n <- nrow(x)
   p <- ncol(x)
   lp <- numeric(n)
   beta <- numeric(p)
   s <- numeric(p)

   for(j in 1:p)
   {
      grad <- sum(x[, j] * (d1phi(lp) - y))
      d2 <- sum(x[,j]^2 * d2phi(lp))
      
      if(grad != 0 && d2 != 0)
      {
	 s[j] <- -grad / d2
	 cat(grad, d2, s[j], "\n")
      	 if(j == 1) {
	    lp <- x[, j] * s[j]
	    cat(">>> lp:", lp[1], "\n");
	 }
      }
   }
   
   #max(abs(s[-1]))
   s
}

(l1max <- getlambda1max(x2, y))
max(abs(l1max[-1]))

stop()

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

