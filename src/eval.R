library(ggplot2)

nexper <- 10
exper <- paste("sim6.", 1:nexper, sep="")
dir <- "~/Software/hapgen_1.3"

n <- 5000
p <- 185805
legend <- sprintf(
   "%s/HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd_legend.txt",
   dir)

#x <- matrix(as.numeric(
#	 readBin(sprintf("%s/%s/sim.bin.t", dir, exper),
#	    what="raw", n=n*(p+1))),
#      nrow=n, byrow=FALSE)
#y <- x[,1]
#x <- x[,-1]


res <- lapply(exper, function(ex) {

   dir <- sprintf("%s/results", ex)

   # Find all files and sort by numerical ordering
   files <- list.files(pattern="^beta_sqrhinge\\.csv\\.", path=dir, full.names=TRUE)
   if(length(files) == 0)
      stop("no files found")
   id <- sapply(strsplit(files, "\\."), function(x) x[4])
   files <- files[order(as.numeric(id))]

   # Get indices of causal SNPs
   snps <- as.character(read.csv(sprintf("%s/loci.txt", ex),
         header=FALSE)[,1])
   snplist <- read.csv(legend, sep="\t")
   pos <- sapply(paste("^", snps, "$", sep=""), grep,
         x=as.character(snplist$position))
   ysnp <- as.numeric((1:p) %in% pos)
   
   b.cd <- sapply(files, function(f) {
      read.csv(f, header=FALSE)[-1,1]
   })
   
   df.cd <- apply(b.cd, 2, function(x) sum(x != 0))

   #res.glm <- lapply(1:ncol(x), function(j) {
   #   cat(j, "\r")
   #   coef(summary(glm(y ~ x[,j], family=binomial)))#[2, 4]
   #})
   #b.glm <- -log10(sapply(res.glm, function(q) if(dim(q)[1] > 1) q[2,4] else 1))
   
   runperf <- function(f)
   {
      s <- system(sprintf("~/Software/perf.src/perf -APR -ROC < %s", f), intern=TRUE)
      p <- as.numeric(sapply(strsplit(s, "[[:space:]]+"), function(x) x[2]))
      names(p) <- sapply(strsplit(s, "[[:space:]]+"), function(x) x[1])
      p
   }
   

   mes <- sapply(1:ncol(b.cd), function(i) {
      f <- sprintf("%s/b.cd.%s", dir, i - 1)
      write.table(cbind(ysnp, abs(b.cd[,i])),
   	 col.names=FALSE, row.names=FALSE, sep="\t",
   	 file=f)
      runperf(f)
   })

   #write.table(cbind(ysnp, b.glm), col.names=FALSE, row.names=FALSE,
   #      file="b.glm", sep="\t")
   #auprc.glm <- runperf("b.glm")

   list(measure=mes, df=df.cd, nsim=length(id))
})

#plot(df.cd, auprc.cd, xlab="# active variables", ylab="MAP", log="x", type="b")
m <- data.frame(
   do.call("rbind", lapply(res, function(r) t(do.call("rbind", r))))
)
m$Sim <- factor(rep(1:length(res), sapply(res, function(x) length(x$df))))

#g.apr <- ggplot(m, aes(x=df, y=APR, colour=Sim)) + geom_point() + geom_line()
#g.apr <- g.apr + scale_x_log10()
#
#g.roc <- ggplot(m, aes(x=df, y=ROC, colour=Sim)) + geom_point() + geom_line()
#g.roc <- g.roc + scale_x_log10()

g.roc <- ggplot(m[m$df > 0, ], aes(x=df, y=ROC))
g.roc <- g.roc + geom_point()
b <- 2^sort(unique(round(log2(m$df[m$df > 0]))))
g.roc <- g.roc + scale_x_log2(labels=b, breaks=b) + stat_smooth()

g.apr <- ggplot(m[m$df > 0, ], aes(x=df, y=APR))
g.apr <- g.apr + geom_point() + scale_x_log2() + stat_smooth()


pdf("results.pdf")
print(g.apr)
print(g.roc)
dev.off()

