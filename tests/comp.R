library(glmnet)
library(ROCR)
library(GeneSets)

n <- 1000
p <- 10000
x <- matrix(rnorm(n * p), n, p)
beta <- numeric(p)
k <- 25
beta[sample(p, k)] <- rnorm(k)
betai <- as.numeric(beta != 0)

y <- ifelse(runif(n) <= plogis(x %*% beta), 1, 0)

p.glm <- lapply(1:p, function(j) {
   cat(j, "\r")
   coef(summary(glm(y ~ x[, j], family=binomial)))
})
pval <- sapply(p.glm, function(q) if(dim(q)[1] > 1) q[2,4] else 1)
lpval <- -log10(pval)

g <- glmnet(x, factor(y), family="binomial")

auprc.glmnet <- apply(as.matrix(coef(g)), 2, function(pr) {
   auprc(abs(pr[-1]), betai, k=500)
})

auprc.univ <- auprc(lpval, betai, k=500)

plot(g$df[-1], auprc.glmnet[-1], type="b", log="x")
abline(v=k, lty=2, lwd=3, col=2)
abline(h=auprc.univ, lty=2)

