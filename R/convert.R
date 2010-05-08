
n <- 1000
p <- 1e2
beta <- rnorm(p + 1)
y <- numeric(n)

f <- file("out.bin", "wb")
for(i in 1:n)
{
   cat(i, "\r")
   x <- rnorm(p)
   y[i] <- ifelse(runif(1) <= drop(plogis(c(1, x) %*% beta)), 1, 0)
   writeBin(y[i], f)
   writeBin(x, f)
}
close(f)
cat("\n")

write.table(y, file="y.csv", sep=",", row.names=FALSE,
      col.names=FALSE)

#stop()
#
#x <- scale(x)
#g <- glm(y ~ x, family=binomial)
#coef(g)
#
#GeneSets::auc(predict(g, type="response"), y, 1, 0)

