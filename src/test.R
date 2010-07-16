# Test CD

library(glmnet)

set.seed(4329210)

n <- 100
p <- 10
w <- sample(0:1, p + 1, replace=TRUE, prob=c(0.9, 0.1)) 
beta <- rnorm(p + 1) * w

x <- matrix(sample(c(0, 1, 2), size=n * p, replace=TRUE), n, p)
y <- ifelse(runif(n) <= plogis(cbind(1, x) %*% beta + rnorm(n, 1, 2)), 1, 0)
table(y)

xs <- scale(x)
z <- cbind(y, xs)
writeBin(as.numeric(z), con="x.bin.t")

g.glmnet <- glmnet(xs, factor(y), family="binomial")
b.glmnet <- as.matrix(coef(g.glmnet))

g.glm <- glm(y ~ xs, family=binomial)
b.glm <- coef(g.glm)

cmd <- sprintf("time ./cd_double -model logistic \\
-f x.bin.t -n %s -p %s \\
-epochs 100 \\
-nl1 100 \\
-beta beta_test.csv -v", n, p)
system(cmd)

files <- list.files(pattern="^beta_test\\.csv\\.")
nm <- as.numeric(sapply(strsplit(files, "\\."), function(x) x[3]))
files <- files[order(nm)]

b.cd <- sapply(files, function(f) {
   read.csv(f, header=FALSE)[,1]
})
df.cd <- apply(b.cd[-1,], 2, function(x) sum(x != 0))

source("../tests/optim.R")
#b.cd2 <- cd2(xs, y, logloss, logd1phi, logd2phi, maxiter=1)

lambda.cd <- read.csv("lambda1path.csv", header=FALSE)[,1]
b.cd2 <- sapply(lambda.cd, cd2, x=xs,
   y=y, lossfunc=logloss, d1phi=logd1phi, d2phi=logd2phi,
   maxiter=100)

par(mfrow=c(1, 3))
matplot(t(b.cd), type="b", pch=21, main="CD-C")
matplot(t(b.cd2), type="b", pch=21, main="CD-R")
matplot(t(b.glmnet), type="b", pch=21, main="glmnet")

stop()

#g1 <- glm(y ~ x, family=binomial)
#b.glm <- coef(g1)
#b.glmnet <- as.matrix(coef(g2))
#b.glmnet.z <- apply(b.glmnet[-1,], 2, function(x) sum(x == 0))
#b.glmnet.1 <- b.glmnet[, which.min((b.glmnet.z - sum(w == 0))^2)]
#
##b.cd2 <- cd2(x, y, logloss, logd1phi, logd2phi, maxiter=100)
#
#res <- cbind(beta=beta, glm=b.glm, cd=b.cd.1, glmnet=b.glmnet.1)
#pairs(res, upper.panel=function(x, y, ...) {
#   points(x, y, ...)
#   abline(0, 1)
#})
#cor(res)
#
#matplot(res, type="l", lty=1)
#
#(rmse.beta <- apply(res, 2, function(x) sqrt(mean(res[,1] - x)^2)))
#(rmse.glm <- apply(res, 2, function(x) sqrt(mean(res[,2] - x)^2)))
#
