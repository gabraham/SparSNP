library(glmnet)
library(GeneSets)

dir <- "sim5"

legend <- "~/Software/hapgen_1.3/HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd_legend.txt"

b.sgd <- read.csv("~/Code/sgd/src/beta_sgd_sim5.csv", header=FALSE)[,1]
#b.cd <- read.csv("~/Code/sgd/src/beta_cd_sim5.csv", header=FALSE)[,1]
x <- matrix(as.numeric(readBin("sim5/sim.bin", what="raw", n=400*185806)),
      nrow=400, byrow=TRUE)
y <- x[,1]
x <- x[,-1]

g <- glmnet(x, factor(y), family="binomial")
b.glmnet <- as.matrix(coef(g))[,11] 

p <- sapply(1:ncol(x), function(j) {
   cat(j, "\r")
   coef(summary(glm(y ~ x[,j], family=binomial)))#[2, 4]
})
p2 <- sapply(p, function(q) if(dim(q)[1] > 1) q[2,4] else  1)
p3 <- -log10(p2)

# Get indices of causal SNPs
snps <- read.csv(sprintf("%s/loci.txt", dir), header=FALSE)[,1]
snplist <- read.csv(legend, sep="\t")
pos <- sapply(snps, grep, x=snplist$position)


ysnp <- as.numeric((1:ncol(x)) %in% pos)

plot(p3)
points(pos, rep(10, 3), col=2, lwd=2)
lines(abs(b.glmnet[-1] * 10), col=3)
lines(abs(b.sgd[-1]) * 10, col=4)

auc(p3, ysnp, 1, 0)
auc(abs(b.sgd[-1]), ysnp, 1, 0)
auc(abs(b.glmnet[-1]), ysnp, 1, 0)

auprc(p3, ysnp)
auprc(abs(b.sgd[-1]), ysnp)
auprc(abs(b.glmnet[-1]), ysnp)

