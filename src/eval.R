library(glmnet)
library(GeneSets)
library(ROCR)

exper <- "sim6"
dir <- "~/Software/hapgen_1.3"

#b.sgd <- read.csv("~/Code/sgd/src/beta_sgd_sim5.csv", header=FALSE)[,1]
#b.gd <- read.csv("~/Code/sgd/src/beta_gd_sim5.csv", header=FALSE)[,1]
#b.cd <- read.csv(sprintf("~/Code/sgd/src/beta_%s_cd.csv", exper), header=FALSE)[,1]
n <- 500
p <- 185805
legend <- sprintf(
   "%s/HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd_legend.txt",
   dir)
#
x <- matrix(as.numeric(
	 readBin(sprintf("%s/%s/sim.bin", dir, exper),
	    what="raw", n=n*(p+1))),
      nrow=n, byrow=TRUE)
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

b.cd <- sapply(0:13, function(i) {
   read.csv(sprintf("beta.csv.%s", i), header=FALSE)[-1,1]
})

df.cd <- apply(b.cd, 2, function(x) sum(x != 0))

acc.cd <- apply(b.cd, 2, function(x) {
   mean(x != 0 & as.numeric(x != 0) == ysnp)
})

g1 <- glmnet(x, factor(y), family="binomial")
b.glmnet <- as.matrix(coef(g1))[-1, ]

acc.glmnet <- apply(b.glmnet, 2, function(x) {
   mean(x != 0 & as.numeric(x != 0) == ysnp)
})

df.glmnet <- g1$df

matplot(cbind(df.cd, df.glmnet), cbind(acc.cd, acc.glmnet), log="",
      type="b", pch=21:22, lty=1, xlab="Non-zero vars",
      ylab="SNP detection accuracy")

#g2 <- glmnet(x, factor(y), family="binomial", alpha=0.5)
#b.glmnet1 <- as.matrix(coef(g1))[,11]
#b.glmnet2 <- as.matrix(coef(g1))[,length(g$lambda)]

#auc.glmnet1 <- apply(as.matrix(coef(g1)), 2, function(p) {
#   auc(abs(p[-1]), ysnp, 1, 0)
#})
#auprc.glmnet1 <- apply(as.matrix(coef(g1)), 2, function(p) {
#   auprc(abs(p[-1]), ysnp, k=n)
#})
#auprc.glmnet1 <- apply(as.matrix(coef(g1)), 2, function(p) {
#   auprc(abs(p[-1]), ysnp, k=500)
#})
#
#p <- lapply(1:ncol(x), function(j) {
#   cat(j, "\r")
#   coef(summary(glm(y ~ x[,j], family=binomial)))#[2, 4]
#})
#p1 <- sapply(p, function(q) if(dim(q)[1] > 1) q[1,4] + q[2,4] else 0)
#p2 <- sapply(p, function(q) if(dim(q)[1] > 1) q[2,4] else  1)
#p3 <- -log10(p2)
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

