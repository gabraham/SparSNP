library(ggplot2)

nexper <- 4
root <- "sim7"
exper <- paste(root, ".", 1:nexper, sep="")
dir <- "~/Software/hapgen_1.3"

n <- 10000
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


# Get indices of causal SNPs
cat("Getting true SNP loci ... ")
snplist <- read.csv(legend, sep="\t")
ysnp <- lapply(exper, function(ex) {
   cat(ex, "\n")
   snps <- as.character(read.csv(sprintf("%s/loci.txt", ex),
         header=FALSE)[,1])
   pos <- sapply(paste("^", snps, "$", sep=""), grep,
         x=as.character(snplist$position))
   as.numeric((1:p) %in% pos)
})
cat("done\n")

runperf <- function(f)
{
   capture.output({s <- system(
	 sprintf("~/Software/perf.src/perf -APR -ROC < %s", f),
	 intern=TRUE)})
   p <- as.numeric(sapply(strsplit(s, "[[:space:]]+"), function(x) x[2]))
   names(p) <- sapply(strsplit(s, "[[:space:]]+"), function(x) x[1])
   p
}

# Analyse CD results
res.cd <- lapply(seq(along=exper), function(k) {
   ex <- exper[[k]]
   dir <- sprintf("%s/results", ex)

   # Find all files and sort by numerical ordering
   files <- list.files(pattern="^beta_sqrhinge\\.csv\\.",
	 path=dir, full.names=TRUE)
   if(length(files) == 0)
      stop("no files found")
   id <- sapply(strsplit(files, "\\."), function(x) x[4])
   files <- files[order(as.numeric(id))]

   b.cd <- sapply(files, function(f) {
      read.csv(f, header=FALSE)[-1,1]
   })
   
   df.cd <- apply(b.cd, 2, function(x) sum(x != 0))

   mes <- sapply(1:ncol(b.cd), function(i) {
      f <- sprintf("%s/b.cd.%s", dir, i - 1)
      if(!file.exists(f)) 
	 write.table(cbind(ysnp[[k]], abs(sign(b.cd[,i]))),
	    col.names=FALSE, row.names=FALSE, sep="\t",
	    file=f)
      cat(f, "\n")
      runperf(f)
   })

   list(measure=mes, df=df.cd, nsim=length(id))
})

m.cd <- data.frame(
   do.call("rbind", lapply(res.cd, function(r) t(do.call("rbind", r))))
)
m.cd$Sim <- factor(
   rep(1:length(res.cd), sapply(res.cd, function(x) length(x$df)))
)

stop()

# Analyse plink results
res.pl <- lapply(seq(along=exper), function(k) {
   ex <- exper[[k]]
   dir <- sprintf("%s/results", ex)
   d <- read.csv(sprintf("%s/plink.assoc.logistic", dir), sep="")
    
   stats <- cbind(coef=d$STAT, logpval=-log10(d$P))

   mes <- sapply(1:ncol(stats), function(i) {
      f <- sprintf("%s/b.cd.plink.%s", dir, colnames(stats)[i])
      if(!file.exists(f))
	 write.table(cbind(ysnp[[k]], abs(stats[,i])),
	    col.names=FALSE, row.names=FALSE, sep="\t",
	    file=f)
      cat(f, "\n")
      runperf(f)
   })
   colnames(mes) <- colnames(stats)
   list(measure=mes)
})

# Both -log10(pval) and STAT yield same AROC/APRC, so take one
m.pl <- data.frame(t(sapply(res.pl, function(x) x[[1]][,1])))
m.pl$df <- 5000 # fake DF, to make the point plot nicely
m.pl$Sim <- 1
m.pl$nsim <- 1
m.pl$Method <- "logistic"

m.cd$Method <- "lasso"

m.comb <- rbind(m.cd, m.pl)
m.comb$Method <- factor(m.comb$Method)

m.comb <- m.comb[m.comb$df > 0 & m.comb$df < 256, ]

maxdf <- round(max(m.comb$df / 10))
m.comb$df_bin <- cut(m.comb$df, breaks=(0:maxdf) * 10,
      labels=(1:maxdf - 1) * 10 + 5)

b <- 2^sort(unique(round(log2(m.comb$df[m.comb$df > 0]))))
g <- ggplot(m.comb, aes(x=df_bin, y=APR)) + geom_boxplot() + geom_point()
#gg <- lapply(c("APR", "ROC"), function(nm) {
#   g <- ggplot(m.comb, aes_string(x="df", y=nm, shape="Method"))
#   g <- g + scale_x_log2(labels=b, breaks=b)
#   g <- g + stat_smooth(method="loess")
#   g <- g + stat_summary(data=subset(m.comb, Method=="logistic"),
#         fun.data=mean_cl_normal, geom="pointrange", colour="blue")
#   g <- g + stat_summary(data=subset(m.comb, Method=="logistic"),
#         fun.data=mean_cl_normal, geom="errorbar", colour="blue")
#   g + geom_point(colour="black")
#})

pdf(sprintf("results_%s.pdf", root))
#for(g in gg) print(g)
print(g)
dev.off()

save.image(file=sprintf("eval_%s.RData", root))

