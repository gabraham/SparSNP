# Test CD

#library(glmnet)

#set.seed(16937)

run <- function(n, p)
{
   y <- 1

   while(length(unique(y)) < 2)
   {
      cat("trying...\n")
      w <- sample(0:1, p + 1, replace=TRUE, prob=c(0.9, 0.1)) 
      beta <- rnorm(p + 1, 0, 0.5) * w
      
      x <- matrix(sample(c(0, 1, 2), size=n * p, replace=TRUE), n, p)
      y <- ifelse(
	 runif(n) <= plogis(cbind(1, x) %*% beta + rnorm(n, 1, 2)), 1, 0)
   }
   
   xs <- scale(x)
   rm(x)
   gc()
   z <- cbind(y, xs)
   writeBin(as.numeric(z), con="x.bin.t")
   
   #g.glmnet <- glmnet(xs, factor(y), family="binomial")
   #b.glmnet <- as.matrix(coef(g.glmnet))
   
   #g.glm <- glm(y ~ xs, family=binomial)
   #b.glm <- coef(g.glm)
   
   cmd <- sprintf("./cd_double -model logistic \\
   -f x.bin.t -n %s -p %s \\
   -epochs 1000 \\
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
   
   lambda.cd <- read.csv("lambda1path.csv", header=FALSE)[,1]
   
   #par(mfrow=c(1, 2))
   #ylim <- range(b.cd, b.glmnet)
   #matplot(df.cd, t(b.cd), type="l", pch=21, main="CD-C", lty=1, ylim=ylim)
   #matplot(g.glmnet$df, t(b.glmnet), type="l", pch=21, main="glmnet",
   #      lty=1, ylim=ylim)
   
   
   runperf <- function(f)
   {
      s <- system(
         sprintf("~/Software/perf.src/perf -APR -ROC < %s", f), intern=TRUE)
      p <- as.numeric(sapply(strsplit(s, "[[:space:]]+"), function(x) x[2]))
      names(p) <- sapply(strsplit(s, "[[:space:]]+"), function(x) x[1])
      p
   }
   
   beta.z <- as.numeric(beta != 0)
   table(beta.z)
   
   mes.cd <- sapply(1:ncol(b.cd), function(i) {
       f <- sprintf("b.cd.%s", i - 1)
       write.table(cbind(beta.z, abs(b.cd[,i])),
    	 col.names=FALSE, row.names=FALSE, sep="\t",
    	 file=f)
       runperf(f)
   })
   
   #mes.glmnet <- sapply(1:ncol(b.glmnet), function(i) {
   #    f <- sprintf("b.cd.%s", i - 1)
   #    write.table(cbind(beta.z, abs(b.glmnet[,i])),
   # 	 col.names=FALSE, row.names=FALSE, sep="\t",
   # 	 file=f)
   #    runperf(f)
   #})
   
   accuracy <- function(a, b)
   {
      mean(a == 0 & b == 0)
   }
   
   acc.cd <- apply(b.cd, 2, accuracy, b=beta.z)
   #acc.glmnet <- apply(b.glmnet, 2, accuracy, b=beta.z)
   
   par(mfrow=c(1, 2))
   plot(df.cd, acc.cd, xlab="# non-zero variables",
      ylab="Accuracy", main="Accuracy")
   matplot(df.cd, t(mes.cd), type="b", pch=21,
      xlab="# non-zero variables", ylab="", main="AUC,APRC")
   par(usr=c(0, 1, 0, 1))
   legend(0.4, 0.5, legend=rownames(mes.cd), col=1:2, lwd=3)
}

rmoutput <- function()
{
   unlink(list.files(pattern="^beta_test\\.csv\\.[[:digit:]]+"))
   unlink(list.files(pattern="^b\\.cd\\.[[:digit:]]+"))
}

rmoutput()
pdf("results1.pdf")
run(n=1e4, p=1e2)
dev.off()

rmoutput()
pdf("results2.pdf")
run(n=1e4, p=1e3)
dev.off()

