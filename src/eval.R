
exper <- "sim6"
dir <- "~/Software/hapgen_1.3"
nruns <- 20

n <- 500
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
snps <- as.character(read.csv(sprintf("%s/%s/loci.txt", dir, exper),
      header=FALSE)[,1])
snplist <- read.csv(legend, sep="\t")
pos <- sapply(paste("^", snps, "$", sep=""), grep,
      x=as.character(snplist$position))
ysnp <- as.numeric((1:p) %in% pos)

b.cd <- sapply(1:nruns - 1, function(i) {
   read.csv(sprintf("beta.csv.%s", i), header=FALSE)[-1,1]
})

df.cd <- apply(b.cd, 2, function(x) sum(x != 0))

#res.glm <- lapply(1:ncol(x), function(j) {
#   cat(j, "\r")
#   coef(summary(glm(y ~ x[,j], family=binomial)))#[2, 4]
#})
#b.glm <- -log10(sapply(res.glm, function(q) if(dim(q)[1] > 1) q[2,4] else 1))

runperf <- function(f)
{
   s <- system(sprintf("~/Software/perf.src/perf -APR < %s", f), intern=TRUE)
   as.numeric(gsub("APR[[:space:]]+", "", s))
}

auprc.cd <- sapply(1:ncol(b.cd), function(i) {
   f <- sprintf("b.cd.%s", i - 1)
   write.table(cbind(ysnp, abs(b.cd[,i])),
	 col.names=FALSE, row.names=FALSE, sep="\t",
	 file=f)
   runperf(f)
})

write.table(cbind(ysnp, b.glm), col.names=FALSE, row.names=FALSE,
      file="b.glm", sep="\t")
auprc.glm <- runperf("b.glm")

plot(df.cd, auprc.cd, xlab="# active variables", ylab="MAP", log="x", type="b")

