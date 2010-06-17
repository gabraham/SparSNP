library(glmnet)

dir <- "sim5"

legend <- "~/Software/hapgen_1.3/HapMap/genotypes_chr1_JPT+CHB_r22_nr.b36_fwd_legend.txt"

b.sgd <- read.csv("~/Code/sgd/src/beta_sgd_sim5.csv", header=FALSE)[,1]
#b.cd <- read.csv("~/Code/sgd/src/beta_cd_sim5.csv", header=FALSE)[,1]
x <- matrix(as.numeric(readBin("sim5/sim.bin", what="raw", n=400*185806)),
      nrow=400, byrow=TRUE)
y <- x[,1]
x <- x[,-1]

g <- glmnet(x, factor(y), family="binomial")

p <- apply(x, 2, function(x) coef(summary(glm(y ~ x, family=binomial)))[2, 4])
p2 <- -log10(p)

# Get indices of causal SNPs
snps <- read.csv(sprintf("%s/loci.txt"), header=FALSE)[,1]
snplist <- read.csv(legend, sep="\t")



