# Process prediction results

library(ggplot2)
library(glmnet)

folds <- read.csv("folds.txt", header=FALSE)[,1]
folds <- unique(sort(folds))
#folds <- 0

sqrhingeloss <- function(p, y)
{
   mean(pmax(0, 1 - p * y)^2)
}

cv.auc <- lapply(folds, function(fold) {
   y <- read.csv(sprintf("y.%02d", fold), header=FALSE)[,1]
   #y <- read.csv("../sim.y", header=FALSE)[,1]
   #y <- (y + 1) / 2

   files <- sort(
      list.files(pattern=sprintf("beta.csv.[[:digit:]]+.%02d.pred", fold))
   )

   pr <- sapply(files, function(f) {
      read.csv(f, header=FALSE)[,1]
   })

   lambda <- read.csv(sprintf("lambda1path.csv.%02d", fold), header=FALSE)[,1]
   nz <- read.csv(sprintf("nonzero.csv.%02d", fold), header=FALSE)[,1]

   # glmnet::auc expects y = 0/1
   res <- apply(pr, 2, auc,
	 y=as.numeric(as.character(factor(y, labels=c(0, 1)))))

   loss <- apply(pr, 2, sqrhingeloss, y=y)

   cbind(lambda=lambda[1:length(res)],
      AUC=res,
      NonZero=nz[1:length(res)],
      loss=loss
   )
})

cv.auc.d <- data.frame(do.call(rbind, cv.auc))

g1 <- ggplot(cv.auc.d, aes(x=NonZero, y=AUC)) 
g1 <- g1 + geom_point() #+ scale_x_log10()
g1 <- g1 + geom_smooth() 

g2 <- ggplot(cv.auc.d, aes(x=lambda, y=AUC)) 
g2 <- g2 + geom_point() #+ scale_x_log10()
g2 <- g2 + geom_smooth() 

g3 <- ggplot(cv.auc.d, aes(x=loss, y=AUC)) 
g3 <- g3 + geom_point() #+ scale_x_log10()
g3 <- g3 + geom_smooth() 

pdf("cv.pdf")
print(g1)
print(g2)
print(g3)
dev.off()


