x <- scan("~/Software/hapgen_1.3/sim4.3/sim.all.g")
dim(x) <- c(1000, 200)
x <- t(x)
y <- scan("~/Software/hapgen_1.3/sim4.3/sim.y")

n <- length(y)
x0 <- numeric(n) + 1


n <- length(y)
p <- ncol(x)
z <- numeric(p)

for(j in 1:p)
{
   x01 <- cbind(1, x[, j])
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
         sum(x01[,1]^2 * w),
         sum(x01[,1] * x01[,2] * w),
         sum(x01[,1] * x01[,2] * w),
         sum(x01[,2]^2 * w)), 2, 2)
      Hinv <- ginv(H)
      s <- Hinv %*% crossprod(x01, lpinv - y) 
      beta <- beta - s
      cat(beta, "\n")
      iter <- iter + 1
   }

   z[j] <- beta[2] / sqrt(diag(Hinv)[2])

}
