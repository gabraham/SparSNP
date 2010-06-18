library(glmnet)
library(GeneSets)

dir <- "sim5"

legend <- "~/Software/hapgen_1.3/HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd_legend.txt.50000"

b.sgd <- read.csv("~/Code/sgd/src/beta_sgd_sim5.csv", header=FALSE)[,1]
b.gd <- read.csv("~/Code/sgd/src/beta_gd_sim5.csv", header=FALSE)[,1]
#b.cd <- read.csv("~/Code/sgd/src/beta_cd_sim5.csv", header=FALSE)[,1]
n <- 1000
p <- 5e4
x <- matrix(as.numeric(readBin("sim5/sim.bin", what="raw", n=n*(p+1))),
      nrow=n, byrow=TRUE)
y <- x[,1]
x <- x[,-1]

# Get indices of causal SNPs
snps <- as.character(read.csv(sprintf("%s/loci.txt", dir), header=FALSE)[,1])
snplist <- read.csv(legend, sep="\t")
pos <- sapply(snps, grep, x=as.character(snplist$position))
ysnp <- as.numeric((1:ncol(x)) %in% pos)


g <- glmnet(x, factor(y), family="binomial")
b.glmnet1 <- as.matrix(coef(g))[,11]
b.glmnet2 <- as.matrix(coef(g))[,length(g$lambda)]

auc.glmnet <- apply(as.matrix(coef(g)), 2, function(p) {
   auc(abs(p[-1]), ysnp, 1, 0)
})
auprc.glmnet <- apply(as.matrix(coef(g)), 2, function(p) {
   auprc(abs(p[-1]), ysnp)
})

p <- lapply(1:ncol(x), function(j) {
   cat(j, "\r")
   coef(summary(glm(y ~ x[,j], family=binomial)))#[2, 4]
})
p1 <- sapply(p, function(q) if(dim(q)[1] > 1) q[1,4] + q[2,4] else 0)
p2 <- sapply(p, function(q) if(dim(q)[1] > 1) q[2,4] else  1)
p3 <- -log10(p2)




res <- cbind(p2, abs(b.glmnet1[-1]), abs(b.glmnet2[-1]), abs(b.sgd[-1]))
matplot(res, type="l", lty=1)
points(pos, rep(max(p3), length(pos)), col=3, lwd=2)

auc(p3, ysnp, 1, 0)
auc(abs(b.gd[-1]), ysnp, 1, 0)
auc(abs(b.sgd[-1]), ysnp, 1, 0)
auc(abs(b.glmnet1[-1]), ysnp, 1, 0)
auc(abs(b.glmnet2[-1]), ysnp, 1, 0)

auprc(p3, ysnp)
auprc(abs(b.gd[-1]), ysnp)
auprc(abs(b.sgd[-1]), ysnp)
auprc(abs(b.glmnet1[-1]), ysnp)
auprc(abs(b.glmnet2[-1]), ysnp)

