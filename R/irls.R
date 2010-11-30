x <- scan("~/Software/hapgen_1.3/sim4.3/sim.all.g")
dim(x) <- c(1000, 200)
x <- t(x)
y <- scan("~/Software/hapgen_1.3/sim4.3/sim.y")

n <- length(y)
x0 <- numeric(n) + 1
x1 <- x[,1]

x01 <- cbind(1, x1)

n <- length(y)
p <- ncol(x)
lp <- numeric(n)
lpinv <- plogis(lp)
beta <- numeric(2)

for(iter in 1:10)
{
   #s0 <- sum(x0 * (lpinv - y)) / sum(x0 ^ 2 * lpinv * (1 - lpinv))
   #s1 <- sum(x1 * (lpinv - y)) / sum(x1 ^ 2 * lpinv * (1 - lpinv))
   #beta <- beta - c(s0, s1)
   #lp <- lp - s0 - x1 * s1
   #lpinv <- plogis(lp)

   lp <- x01 %*% beta
   lpinv <- plogis(lp)
   w <- lpinv * (1 - lpinv)
   H <- matrix(c(
      sum(x0^2 * w),
      sum(x0 * x1 * w),
      sum(x1 * x0 * w),
      sum(x1^2 * w)), 2, 2)
   beta <- beta - solve(H) %*% crossprod(x01, lpinv - y) 
   cat(beta, "\n")
}


