# Test CD

library(glmnet)

#set.seed(43249210)

n <- 1000
p <- 100
w <- sample(0:1, p + 1, replace=TRUE, prob=c(0.9, 0.1)) 
beta <- rnorm(p + 1) * w

x <- matrix(sample(c(0, 1, 2), size=n * p, replace=TRUE), n, p)
y <- ifelse(runif(n) <= plogis(cbind(1, x) %*% beta + rnorm(n, 1, 2)), 1, 0)
table(y)

z <- cbind(y, x)

writeBin(as.raw(z), con="x.bin.t")

nl1 <- 50
cmd <- sprintf("time ./cd -model logistic -f x.bin.t -n %s -p %s -beta beta.csv \\
-epochs 100 -thresh 1e-4 -v -nl1 %s", n, p, nl1)
system(cmd)

b.cd <- sapply(1:nl1 - 1, function(i) {
   read.csv(sprintf("beta.csv.%s", i), header=FALSE)[,1]
})
b.cd.z <- apply(b.cd[-1,], 2, function(x) sum(x == 0))
g1 <- glm(y ~ x, family=binomial)
g2 <- glmnet(x, factor(y), family="binomial", lambda=0.2)
b.glm <- coef(g1)
b.glmnet <- as.numeric(coef(g2))

source("../tests/optim.R")
b.cd2 <- cd2(x, y, logloss, logd1phi, logd2phi, maxiter=100)

res <- cbind(beta=beta, glm=b.glm, cd=b.cd, cd2=b.cd2, glmnet=b.glmnet)
pairs(res, upper.panel=function(x, y, ...) {
   points(x, y, ...)
   abline(0, 1)
})
cor(res)

matplot(res, type="l", lty=1)

(rmse.beta <- apply(res, 2, function(x) sqrt(mean(res[,1] - x)^2)))
(rmse.glm <- apply(res, 2, function(x) sqrt(mean(res[,2] - x)^2)))

