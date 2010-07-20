# Test CD

library(ggplot2)

#library(glmnet)

#set.seed(16937)

runperf <- function(f)
{
   s <- system(
      sprintf("~/Software/perf.src/perf -APR -ROC < %s", f), intern=TRUE)
   p <- as.numeric(sapply(strsplit(s, "[[:space:]]+"), function(x) x[2]))
   names(p) <- sapply(strsplit(s, "[[:space:]]+"), function(x) x[1])
   p
}

run <- function(n, p, nsim=50)
{
   sim <- function()
   {
      y <- 1

      while(length(unique(y)) < 2)
      {
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
      -epochs 100 \\
      -nl1 100 -v \\
      -beta beta_test.csv 2>&1 > /dev/null", n, p)
      out <- capture.output(system(cmd))
      
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

      cbind(df=df.cd, acc=acc.cd, mes=t(mes.cd))
   }

   res <- lapply(1:nsim, function(i) {
      cat(i, "")
      sim()
   })

   res2 <- data.frame(do.call("rbind", res), check.names=FALSE,
	 check.rows=FALSE)
   res2$Sim <- rep(1:nsim, sapply(res, nrow))
   colnames(res2) <- c("DF", "ACC", "APRC", "AROC", "Sim")
   #res2$DF <- factor(res2$DF)

   res2
}
   

rmoutput <- function()
{
   unlink(list.files(pattern="^beta_test\\.csv\\.[[:digit:]]+"))
   unlink(list.files(pattern="^b\\.cd\\.[[:digit:]]+"))
}

do_plot <- function(r)
{
   g1 <- ggplot(r, aes(x=DF, y=AROC)) 
   g1 <- g1 + ylim(0, 1)
   g1 <- g1 + stat_summary(fun.y="mean_cl_normal",
      geom="errorbar", fun.ymin=min, fun.ymax=max)
   g1 <- g1 + stat_summary(fun.y=mean, geom="point")

   g2 <- ggplot(r, aes(x=DF, y=APRC)) 
   g2 <- g2 + ylim(0, 1)
   g2 <- g2 + stat_summary(fun.y="mean_cl_normal",
      geom="errorbar", fun.ymin=min, fun.ymax=max)
   g2 <- g2 + stat_summary(fun.y=mean, geom="point")

   grid.newpage()
   pushViewport(viewport(layout=grid.layout(1, 2)))

   print(g1, vp=viewport(layout.pos.row=1, layout.pos.col=1))
   print(g2, vp=viewport(layout.pos.row=1, layout.pos.col=2))
}

res <- lapply(4:10, function(i) {
   rmoutput()
   cat(i, "... ")
   r <- run(n=1000, p=2^i, nsim=50)
   cat("\n")
   r
})

pdf("results_sim.pdf", width=10)
for(r in res) do_plot(r)
dev.off()

