
n <- 1000
p <- 1e4

x <- matrix(rnorm(n * p), n, p)
y <- ifelse(runif(n) <= plogis(rowSums(cbind(1, scale(x)))), 1, 0)

d <- as.numeric(t(cbind(y, x)))

f <- file("out.bin", "wb")
writeBin(d, f)
close(f)

#stop()
#
#x <- scale(x)
#g <- glm(y ~ x, family=binomial)
#coef(g)
#
#GeneSets::auc(predict(g, type="response"), y, 1, 0)

