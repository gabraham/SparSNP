# Testing the C sgd code

library(glmnet)
library(penalized)

source("../R/convert.R")

pdf("testsgd.pdf")

# n >> p unpenalised logistic regression
s <- simulate(outfile="sim.bin", n=1000, p=20)
g <- glm(s$y ~ s$x, family=binomial)
system("../src/sgd -m logistic -f sim.bin -n 1000 -p 20 -e 1e2 -v -s 0.5 -t 1e-6")
beta <- read.csv("beta.csv", header=FALSE)[,1]
yhat <- read.csv("pred.csv", header=FALSE)[,1]
#compare(beta1=coef(g), yhat1=predict(g, type="response"),
#   beta2=beta, yhat2=yhat, xlab="glm", ylab="sgd")
d.beta <- cbind(sgd=beta, glm=as.numeric(coef(g)), beta=s$beta)
d.pred <- cbind(sgd=yhat, glm=predict(g, type="response"))
pairs(d.beta, main="Weights")
pairs(d.pred, main="Predictions")

# n << p L2 penalised logistic regression
n <- 1000
p <- 5000
s <- simulate(outfile="sim.bin", n=n, p=p, noise=rnorm(n, 0, 0.2))
lambda2 <- 1e-6

pn <- penalized(s$y, s$x, lambda2=lambda2, model="logistic")

maxit <- c(10, 50, 100, 200)

l1 <- lapply(maxit, function(e) {
   system(sprintf("../src/sgd -m logistic -f sim.bin -n %d -p %d e3 -e %d \\
   -s 1e-6 -t 1e-6 -l2 %.5f -t 1e-9", n, p, e, lambda2))
   list(
      beta=read.csv("beta.csv", header=FALSE)[,1],
      yhat=read.csv("pred.csv", header=FALSE)[,1]
   )
})
l1.beta <- do.call("cbind", lapply(l1, function(x) x$beta))
l1.yhat <- do.call("cbind", lapply(l1, function(x) x$yhat))
colnames(l1.beta) <- paste("sgd", maxit)
colnames(l1.yhat) <- paste("sgd", maxit)

cat("glmnet...\n")
l2 <- lapply(maxit, function(e) {
   g <- glmnet(s$x, factor(s$y), family="binomial",
	 alpha=0, lambda=lambda2, maxit=e)
   list(
      beta=as.numeric(coef(g)),
      yhat=drop(predict(g, s$x, type="response"))
   )
})
l2.beta <- do.call("cbind", lapply(l2, function(x) x$beta))
l2.yhat <- do.call("cbind", lapply(l2, function(x) x$yhat))
colnames(l2.beta) <- paste("glmnet", maxit)
colnames(l2.yhat) <- paste("glmnet", maxit)

d.beta <- cbind(l1.beta, l2.beta, penalized=as.numeric(coef(pn, "all")),
      beta=s$beta)
d.pred <- cbind(l1.yhat, l2.yhat, penalized=as.numeric(predict(pn, s$x)))

pairs(d.beta, main="Weights")
pairs(d.pred, main="Predictions")

dev.off()

