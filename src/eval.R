library(glmnet)
library(GeneSets)
library(ROCR)

exper <- "sim6"
dir <- "~/Software/hapgen_1.3"

#b.sgd <- read.csv("~/Code/sgd/src/beta_sgd_sim5.csv", header=FALSE)[,1]
#b.gd <- read.csv("~/Code/sgd/src/beta_gd_sim5.csv", header=FALSE)[,1]
#b.cd <- read.csv(sprintf("~/Code/sgd/src/beta_%s_cd.csv", exper), header=FALSE)[,1]
n <- 1000
p <- 185805
legend <- sprintf(
   "%s/HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd_legend.txt",
   dir)
#
x <- matrix(as.numeric(
	 readBin(sprintf("%s/%s/sim.bin.t", dir, exper),
	    what="raw", n=n*(p+1))),
      nrow=n, byrow=FALSE)
y <- x[,1]
x <- x[,-1]

################################################################################
# Compare coordinate descent with glm() on small data p << N

#logloss <- function(x, y, beta)
#{
#   lp <- x %*% beta
#   log(1 + exp(lp)) - y * lp
#}
#
#system.time({
#   g1 <- glm(y ~ x, family=binomial())
#})
#
#dir <- sprintf("~/Software/hapgen_1.3/%s", exper)
#system(sprintf("~/Code/cd/src/cd -model logistic -f \\
#%s/sim.bin.t -n %s -p %s -epochs 100 \\
#-v -beta %s/beta_%s_cd.csv", dir, n, p, dir, exper))
#b.cd <- read.csv(
#   sprintf("%s/beta_%s_cd.csv", dir, exper), header=FALSE)[,1]
#
#
#
#stop()


# Get indices of causal SNPs
snps <- as.character(read.csv(sprintf("%s/%s/loci.txt", dir, exper),
      header=FALSE)[,1])
snplist <- read.csv(legend, sep="\t")
pos <- sapply(paste("^", snps, "$", sep=""), grep,
      x=as.character(snplist$position))
ysnp <- as.numeric((1:ncol(x)) %in% pos)

b.cd <- sapply(0:8, function(i) {
   read.csv(sprintf("beta.csv.%s", i), header=FALSE)[-1,1]
})

df.cd <- apply(b.cd, 2, function(x) sum(x != 0))

acc.cd <- apply(b.cd, 2, function(x) {
   mean(x != 0 & as.numeric(x != 0) == ysnp)
})

cd.auprc <- sapply(1:ncol(b.cd), function(i) auprc(abs(b.cd[,i]), ysnp))

res.glm <- lapply(1:ncol(x), function(j) {
   cat(j, "\r")
   coef(summary(glm(y ~ x[,j], family=binomial)))#[2, 4]
})
b.glm <- -log10(sapply(res.glm, function(q) if(dim(q)[1] > 1) q[2,4] else 1))
b.auprc <- auprc(b.glm, ysnp)
#



#res <- cbind(p2, abs(b.glmnet1[-1]), abs(b.glmnet2[-1]), abs(b.sgd[-1]))
#matplot(res, type="l", lty=1)
#points(pos, rep(max(p3), length(pos)), col=3, lwd=2)

#auc(p3, ysnp, 1, 0)
#auc(abs(b.cd[-1]), ysnp, 1, 0)
#auc(abs(b.gd[-1]), ysnp, 1, 0)
#auc(abs(b.sgd[-1]), ysnp, 1, 0)
#auc(abs(b.glmnet1[-1]), ysnp, 1, 0)
#auc(abs(b.glmnet2[-1]), ysnp, 1, 0)

#auprc(p3, ysnp, k=500)
#auprc(abs(b.cd[-1]), ysnp, k=500)

#auprc(abs(b.gd[-1]), ysnp, k=500)
#auprc(abs(b.sgd[-1]), ysnp, k=500)
#auprc(abs(b.glmnet1[-1]), ysnp, k=500)
#auprc(abs(b.glmnet2[-1]), ysnp, k=500)
#
#
#par(mfrow=c(2, 3))
#l <- list(GD=b.gd[-1], SGD=b.sgd[-1], Univ=p3,
#   glmnet1=b.glmnet1[-1], glmnet2=b.glmnet2[-1])
#for(i in seq(along=l))
#{
#   plot(performance(prediction(l[[i]], ysnp), "prec", "rec"), ylim=c(0, 1),
#      main=names(l)[i])
#}
#

#plot(g1$df[-1], auprc.glmnet1[-1], type="b", log="x",
#      xlab="Nonzero variables", ylab="AUPRC")

