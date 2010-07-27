library(GeneSets)
library(glmnet)
library(ElemStatLearn)

#source("~/Code/sgd/R/gmatrix.R")
#source("~/Code/sgd/R/sgd.R")

#test.multinom <- function()
#{
#   set.seed(10278)
#   K <- 5
#   N <- 1000
#   p <- 10
#   x <- matrix(rnorm(N * p), N, p)
#   #beta <- rnorm(p)
#   y <- sample(1:K, N, TRUE)
#   ##g <- glmnet(x, factor(y), family="multinomial", alpha=0, lambda=0)
#   ##pr <- predict(g, x)
#   #Y <- sapply(1:K, function(i) as.numeric(y == i))
#   #stepsize <- 1e-3
#
#   #B <- matrix(0, K, p)
#   #loss <- numeric(N)
#   #for(i in 1:N)
#   #{
#   #   loss[i] <- mlogloss(x[i, , drop=FALSE], Y[i,,drop=FALSE], B) 
#   #   cat("loss:", loss[i], "\n")
#   #   B <- B - stepsize * mlogdloss(x[i,,drop=FALSE], Y[i,,drop=FALSE], B) 
#   #}
# 
#   #q <- exp(tcrossprod(x, B))
#   #p <- q / rowSums(q)
#   #mean(apply(p, 1, which.max) != y)
#   #mean(apply(pr, 1, which.max) != y)
#   fname <- tempfile()
#   out <- file(fname, "wb")
#   writeBin(as.numeric(t(x)), out)
#   close(out)
#
#   m1 <- sgd(fname, y=y, p=p, model="multinomial", maxepochs=300)
#   m2 <- glmnet(x, factor(y), family="multinomial", alpha=0, lambda=0)
#   pr <- exp(tcrossprod(cbind(1, x), m1)) 
#   p1 <- pr / rowSums(pr)
#   p2 <- predict(m2, x, type="response")[,,1]
#   list(sgd=m1, glmnet=m2, x=x, y=y, p1=p1, p2=p2)
#}


#set.seed(643694)
#n <- 500L
#p <- 10L
#x <- matrix(rnorm(n * p), n, p)
#beta <- cbind(c(1, rnorm(p, 0, 1) * sample(0:1, p, replace=TRUE)))
#y <- cbind(1, x) %*% beta + rnorm(n)

library(ElemStatLearn)
y <- prostate$lpsa[prostate$train]
x <- base::scale(as.matrix(prostate[, 1:8]))[prostate$train, ]
n <- nrow(x)
p <- ncol(x)

g1 <- gmatrixMem(x, nrow=n, ncol=p, companions=list(y=cbind(y)))

pen <- 2^(9:-9)

l1 <- sapply(pen, function(l) {
   reset(g1)
   cat("lambda1:", l, "\n")
   m <- sgd.gmatrix(g=g1, model="linear", maxepochs=50, stepsize=5e-4,
      threshold=1e-6, lambda1=l, verbose=TRUE)
   m@B
})

l2 <- sapply(pen, function(l) {
   reset(g1)
   cat("lambda2:", l, "\n")
   m <- sgd.gmatrix(g=g1, model="linear", maxepochs=1e3, stepsize=5e-4,
      threshold=1e-9, lambda2=l, verbose=TRUE)
   m@B
})


rownames(l1) <- rownames(l2) <- c("(Intercept)", colnames(x))

#(cf <- cbind(true=beta, l2=m1@B, l1=l1))
#par(mfrow=c(1, 3))
#matplot(cf, type="b", lty=1)
#abline(h=0, lty=2)
#par(usr=c(0, 1, 0, 1))
#legend(0.8, 0.7, legend=c("True", "L2", "L1"), lwd=3, col=1:3)

l <- coef(lm(y ~ x))

pdf("prostate.pdf", width=12)
par(mfrow=c(1, 2))

#l1[l1 == 0] <- as.numeric(NA)
ylim <- range(l1[-1,], l2[-1,], l[-1], na.rm=TRUE)


matplot(1:length(pen), t(l1)[,-1], type="b", main="L1",
      xlab="Penalisation", ylab="Coef", pch=20, lty=1,
      ylim=ylim, xlim=c(1, length(pen) + 5), col=1:p, cex=0.95)
points(rep(length(pen)+1, p), l[-1], col=1:p, pch=22, cex=0.95)
abline(h=0, lty=3)
text(length(pen) + 3, l[-1], labels=rownames(l1)[-1], col=1:p)

matplot(1:length(pen), t(l2)[,-1], type="b", main="L2",
      xlab="Penalisation", ylab="Coef", pch=20, lty=1,
      ylim=ylim, xlim=c(1, length(pen) + 5), col=1:p, cex=0.95)
points(rep(length(pen)+1, p), l[-1], col=1:p, pch=22, cex=0.95)
abline(h=0, lty=3)
text(length(pen) + 3, l[-1], labels=rownames(l1)[-1], col=1:p)

dev.off()

stop()

# Classification
n <- 100L
p <- 10L
x <- matrix(rnorm(n * p), n, p)
beta <- rnorm(p + 1)
g1 <- gmatrixMem(x, nrow=n, ncol=p, companions=list(y=cbind(y)))
y <- ifelse(runif(nrow(x)) < plogis(cbind(1, x) %*% beta), 1, -1)
y2 <- ifelse(y == -1, 0, 1)

g <- list(g1, g1, g1)
g[[1]]@companions$y <- cbind(y2)
g[[2]]@companions$y <- cbind(y)
g[[3]]@companions$y <- cbind(y)

models <- c("logistic", "boosting", "hinge") 

m <- lapply(seq(along=models), function(i) {
   reset(g[[i]])
   sgd.gmatrix(g[[i]], model=models[[i]], stepsize=1e-3,
	 maxepochs=1e3, threshold=1e-8)@B
})

err <- sapply(seq(along=m), function(i) {
   reset(g[[i]])
   mean(sign(predict(new("sgd", B=cbind(m[[i]])), g[[i]])) != y)
})

reset(g1)
cv <- crossval(g1, nfolds=5, nreps=10, nsamples=n,
      model="linear", maxepoch=20, stepsize=1e-2,
      eval=function(x, y) mean((x-y)^2), verbose=FALSE)

stop()

##
##l1 <- lm(y ~ x)
##
##par(mfrow=c(1, 3))
##
##plot(y, type="b")
##pr1 <- cbind(1, x) %*% m1 
##lines(pr1, col=2, type="b")
##lines(fitted(l1), col=3, type="b")
##
##plot(fitted(l1), pr1)
##abline(0, 1)
##
##d <- cbind(beta, coef(l1), m1)
##matplot(d, type="b", lty=1)
#
#y <- ifelse(runif(n) < plogis(cbind(1, x) %*% beta), 1, 0)
##g2 <- gmatrixMem(x, nrow=n, ncol=p, companions=list(y=y))
##m2 <- sgd.gmatrixMem(g=g2, model="logistic", saveloss=TRUE)
#l2 <- glm(y ~ x, family=binomial)
##pr2 <- 1 / (1 + exp(-cbind(1, x) %*% m2))
#f2 <- fitted(l2)
##
##auc(pr2, y, 1, 0)
##auc(fitted(l2), y, 1, 0)
#
###par(mfrow=c(1, 2))
##
##matplot(cbind(beta, m2, coef(l2)), lty=1, type="b", pch=21)
##
##matplot(cbind(
##      cbind(1, x) %*% beta,
##      predict(l2, type="link"),
##      cbind(1, x) %*% m2
##   ),
##   type="b", lty=1, pch=21
##)
#
#Y <- sapply(0:1, function(i) as.numeric(y == i))
#
##reset(g2)
#g2 <- gmatrixMem(x, nrow=n, ncol=p, companions=list(y=y))
#m2 <- sgd.gmatrix(g=g2, model="logistic", saveloss=TRUE, maxepochs=1e2,
#      anneal=0)
#
#cat("#########################\n")
##reset(g2)
#g2 <- gmatrixMem(x, nrow=n, ncol=p, companions=list(y=y))
#m3 <- sgd.gmatrix(g=g2, model="multinomial", saveloss=TRUE, maxepochs=1e2,
#      anneal=0)
#
#pr2 <- 1 / (1 + exp(-cbind(1, x) %*% coef(m2)))
#pr3 <- exp(cbind(1, x) %*% coef(m3)) / rowSums(exp(cbind(1, x) %*% coef(m3)))
#
#n1 <- glmnet(x, factor(y), family="binomial", lambda=0)
#n2 <- glmnet(x, factor(y), family="multinomial", lambda=0)
#pn1 <- glmnet:::predict.lognet(n1, x, type="response")
#
#fname <- "test.dat"
#f <- file(fname, "wb")
#writeBin(as.numeric(t(x)), f)
#close(f)
#
#g4 <- gmatrixDisk(fname, nrow=n, ncol=p, companions=list(y=y))
#m4 <- sgd.gmatrix(g=g4, model="multinomial", saveloss=TRUE, maxepochs=1e2,
#      anneal=0)
#pr4 <- exp(cbind(1, x) %*% coef(m4)) / rowSums(exp(cbind(1, x) %*% coef(m4))) 
#
#
#b1 <- b2 <- cbind(numeric(p + 1))
#X <- cbind(1, x)
#for(k in 1:200)
#{
#   for(i in 1:n)
#   {
#      grad1 <- X[i,] * (exp(t(b1) %*% X[i,]) / (1 + exp(t(b1) %*% X[i,])) - y[i])
#      grad2 <- logdloss(X[i,,drop=FALSE], y[i,,drop=FALSE], b2)
#      b1 <- b1 - 1e-3 * grad1
#      b2 <- b2 - 1e-3 * grad2
#   }
#   #cat(k, b1, ":", b2, "\n")
#}
#
#pra1 <- 1 / (1 + exp(-X %*% b1))
#pra2 <- 1 / (1 + exp(-X %*% b2))
#
#pdf("results.pdf")
#par(mfrow=c(2, 3), pty="s")
#plot(pn1, pr2, xlab="glmnet", ylab="logistic")
#abline(0, 1, lty=2)
#plot(pn1, pr3[, 2], xlab="glmnet", ylab="multinomial")
#abline(0, 1, lty=2)
#plot(f2, pr2, xlab="glm", ylab="logistic")
#abline(0, 1, lty=2)
#plot(f2, pr3[, 2], xlab="glm", ylab="multinomial")
#abline(0, 1, lty=2)
#plot(f2, pra1, xlab="glm", ylab="adhoc")
#abline(0, 1)
#plot(pr2, pra1, xlab="logistic", ylab="adhoc")
#abline(0, 1)
#
#par(mfrow=c(1, 1))
#matplot(cbind(
#      true=beta,
#      glm=coef(l2),
#      sgdlogis=coef(m2),
#      sgdmultinom=coef(m3)[,2],
#      glmnetlogis=as.matrix(coef(n1))[,1],
#      glmnetmultunom=as.matrix(coef(n2)[[2]])[,1]
#   ),
#   lwd=3, cex=3, type="b"
#)
#
#dev.off()
#
#reset(g2)
#ppr2 <- predict(m2, g2, type="response")
#mean((ppr2-pr2)^2)
#
#reset(g2)
#ppr3 <- predict(m3, g2, type="response")
#mean((ppr3-pr3)^2)
#
#reset(g4)
#ppr4 <- predict(m4, g4, type="response")
#mean((ppr4-pr4)^2)

################################################################################
# Multiclass

#x <- as.matrix(vowel.test[, -1])
#g1 <- gmatrixMem(x, nrow=nrow(x), ncol=ncol(x),
#      companions=list(y=cbind(vowel.test$y)))
#m1 <- sgd.gmatrix(g=g1, model="multinomial") 
#reset(g1)
#p1 <- predict(m1, g1, type="response")
#mean(apply(p1, 1, which.max) != vowel.test$y)
#
#m2 <- glmnet(x, factor(vowel.test$y), family="multinomial")
#p2 <- predict(m2, x, type="response")
#mean(apply(p2, 1, which.max) != vowel.test$y)

# From http://www.jstatsoft.org/v33/i01/
# Ramaswamy 14-class cancer dataset
load("~/Data/Ramaswamy/Ramaswamy.RData")
x <- rbind(Ramaswamy$Train$x, Ramaswamy$Test$x)
y <- c(Ramaswamy$Train$y, Ramaswamy$Test$y)
Y <- sapply(sort(unique(y)), function(i) as.numeric(y == i))

#x <- x[1:20, ]
#y <- y[1:20]
#Y <- sapply(sort(unique(y)), function(i) as.numeric(y == i))

g1 <- gmatrixMem(x, nrow=nrow(x), ncol=ncol(x),
      companions=list(y=cbind(y)))
#m1 <- sgd.gmatrix(g=g1, model="multinomial", lambda2=0.03, maxepochs=5,
#      saveloss=TRUE, threshold=1e-6) 
#reset(g1)
#p1 <- predict(m1, g1, type="response")
#(err1 <- mean(apply(p1, 1, which.max) != y))
#mlogloss(cbind(1, x), Y, m1@B)

cv <- crossval(g1, nsamples=nrow(x), nfolds=3, model="multinomial",
   lambda2=0.03, maxepochs=10, verbose=FALSE,
   eval=function(pr, y) {
      mean(apply(pr, 1, which.max) != y)	 
   }
)


nfolds <- 3
nreps <- 1
yf <- factor(y)
nsamples <- nrow(x)
cvglmnet <- lapply(1:nreps, function(rep) {
   folds <- sample(nfolds, size=nsamples, replace=TRUE)
   lapply(1:nfolds, function(fold) {
      m <- glmnet(x[folds != fold, ], yf[folds != fold], family="multinomial")
      pr.train <- predict(m, x[folds != fold, ], type="response")
      pr.test <- predict(m, x[folds == fold, ], type="response")

      res.train <- sapply(1:100, function(i) mean(apply(pr.train[,,i], 1,
	       which.max) != y[folds != fold]))

      res.test <- sapply(1:100, function(i) mean(apply(pr.test[,,i], 1,
	       which.max) != y[folds == fold]))
      #cat("train:", res.train, "test:", res.test, "\n")
      cbind(train=res.train, test=res.test)
   })
})


cvglm <- cv.glmnet(x, factor(y), family="multinomial",
      alpha=0, lambda=0.03, nfolds=3, type="class")

stop()

## Warm start -- we need to supply both the best loss and its coefs and the
# last loss and its coefs so that the fitting can continue from where we
# stopped...
#
#reset(g1)
#m2 <- sgd.gmatrix(g=g1, B=m1@B, loss=m1@loss, model="multinomial",
#      lambda2=0.03, maxepochs=5, stepsize=m1@stepsize)
#reset(g1)
#p2 <- predict(m2, g1, type="response")
#(err2 <- mean(apply(p2, 1, which.max) != y))
#mlogloss(cbind(1, x), Y, m2@B)

# Long run
reset(g1)
m3 <- sgd.gmatrix(g=g1, model="multinomial", lambda2=0.03, maxepochs=20,
      saveloss=TRUE, threshold=1e-6) 
reset(g1)
p3 <- predict(m3, g1, type="response")
(err3 <- mean(apply(p3, 1, which.max) != y))
mlogloss(cbind(1, x), Y, m3@B)



#m3 <- glmnet(x, factor(y), family="multinomial", lambda=0.03, alpha=0)
#p3 <- predict(m3, x, type="response")
#err3 <- sapply(1:dim(p3)[3], function(i) {
#   mean(apply(p3[,,i], 1, which.max) != y)
#})
#mlogloss(cbind(1, x), Y, m2@B)


