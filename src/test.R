# Test CD

#library(glmnet)

#set.seed(16937)

run <- function(n, p, nsim=3)
{
   sim <- function()
   {
      y <- 1

      while(length(unique(y)) < 2)
      {
         cat("trying...\n")
         beta.z <- 0
         while(all(beta.z == 0))
            beta.z <- sample(0:1, p + 1, replace=TRUE, prob=c(0.9, 0.1)) 
         beta <- rnorm(p + 1, 0, 0.5) * beta.z
         
         x <- matrix(sample(c(0, 1, 2), size=n * p, replace=TRUE), n, p)
         y <- ifelse(
            runif(n) <= plogis(cbind(1, x) %*% beta + rnorm(n, 1, 2)), 1, 0)
      }

      xs <- scale(x)
      rm(x)
      gc()
      z <- cbind(y, xs)
      writeBin(as.numeric(z), con="x.bin.t")
      
      cmd <- sprintf("./cd_double -model logistic \\
      -f x.bin.t -n %s -p %s \\
      -epochs 2000 \\
      -nl1 100 -v \\
      -beta beta_test.csv", n, p)
      system(cmd)
      
      files <- list.files(pattern="^beta_test\\.csv\\.")
      if(length(files) == 0)
         return()

      nm <- as.numeric(sapply(strsplit(files, "\\."), function(x) x[3]))
      files <- files[order(nm)]
      
      b.cd <- sapply(files, function(f) {
         read.csv(f, header=FALSE)[,1]
      })
      df.cd <- apply(b.cd[-1,], 2, function(x) sum(x != 0))
      
      lambda.cd <- read.csv("lambda1path.csv", header=FALSE)[,1]
      
      runperf <- function(f)
      {
         s <- system(
            sprintf("~/Software/perf.src/perf -APR -ROC < %s", f), intern=TRUE)
         p <- as.numeric(sapply(strsplit(s, "[[:space:]]+"), function(x) x[2]))
         names(p) <- sapply(strsplit(s, "[[:space:]]+"), function(x) x[1])
         p
      }
      
      mes.cd <- sapply(1:ncol(b.cd), function(i) {
          f <- sprintf("b.cd.%s", i - 1)
          write.table(cbind(beta.z, abs(b.cd[,i])),
       	 col.names=FALSE, row.names=FALSE, sep="\t",
       	 file=f)
          runperf(f)
      })
      
      accuracy <- function(a, b)
      {
         mean(a == 0 & b == 0)
      }
      
      acc.cd <- apply(b.cd, 2, accuracy, b=beta.z)

      list(df=df.cd, acc=acc.cd, mes=mes.cd)
   }

   res <- replicate(nsim, sim(), simplify=FALSE)

   browser()
   
   par(mfrow=c(1, 2))
   plot(df.cd, acc.cd, xlab="# non-zero variables",
      ylab="Accuracy", main="Accuracy", ylim=c(0, 1))
   matplot(df.cd, t(mes.cd), type="b", pch=21,
      xlab="# non-zero variables", ylab="", main="AUC,APRC", ylim=c(0, 1))
   par(usr=c(0, 1, 0, 1))
   legend(0.7, 0.2, legend=rownames(mes.cd), col=1:2, lwd=3)
}

rmoutput <- function()
{
   unlink(list.files(pattern="^beta_test\\.csv\\.[[:digit:]]+"))
   unlink(list.files(pattern="^b\\.cd\\.[[:digit:]]+"))
}

pdf("results_sim.pdf", width=10)

for(i in 2:13)
{
   rmoutput()
   run(n=1000, p=2^i)
}

dev.off()

