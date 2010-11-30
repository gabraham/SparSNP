x <- scan("~/Software/hapgen_1.3/sim4.3/sim.all.g")
dim(x) <- c(1000, 200)
x <- t(x)
y <- scan("~/Software/hapgen_1.3/sim4.3/sim.y")

n <- length(y)
x0 <- numeric(n) + 1
x1 <- x[,3]

x01 <- cbind(1, x1)

n <- length(y)
p <- ncol(x)
lp <- numeric(n)
lpinv <- plogis(lp)
beta <- numeric(2)
iter <- 1
s <- numeric(2) + 10

while(any(abs(s) >= 1e-9) && iter <= 100)
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
   s <- solve(H) %*% crossprod(x01, lpinv - y) 
   beta <- beta - g * s
   cat(beta, "\n")
   iter <- iter + 1
}


