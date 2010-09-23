# Process prediction results

library(ggplot2)
library(glmnet)

folds <- read.csv("folds.txt", header=FALSE)[,1]
folds <- unique(sort(folds))

cv.auc <- sapply(folds, function(fold) {
   y <- read.csv(sprintf("y.%02d", fold), header=FALSE)[,1]
   y01 <- (y + 1) / 2

   files <- sort(
      list.files(pattern=sprintf("beta.csv.[[:digit:]]+.%02d.pred", fold))
   )

   pr <- sapply(files, function(f) {
      read.csv(f, header=FALSE)[,1]
   })

   lambda <- read.csv(sprintf("lambda1path.csv.%02d", fold), header=FALSE)[,1]
   nz <- read.csv(sprintf("nonzero.csv.%02d", fold), header=FALSE)[,1]

   res <- apply(pr, 2, auc, y=y01)
   
   cbind(lambda=lambda[1:length(res)], AUC=res, NonZero=nz[1:length(res)])
})

cv.auc.d <- data.frame(do.call(rbind, cv.auc))

g <- ggplot(cv.auc.d, aes(x=NonZero, y=AUC)) 
g <- g + geom_point() + scale_x_log2()
g <- g + geom_smooth() 

pdf("cv.pdf")
print(g)
dev.off()


