# Genetic variance explained
# Wray et al "Genetic interpretation of Area under the ROC Curve in Genomic
# Profiling", PLoS Genetics, 2:e1000864, 2010
#
#
# K: prevelance > 0
# h2l: heritability
# auc: Area under ROC Curve
# 
# Must specify K and at least one of h2l, auc
#
# Output:
#   K:      prevalence as supplied
#   h2l:    genetic heritability as supplied
#   aucmax: maximum AUC given K and h2l (explains all genetic variance)
#   auc:    AUC as supplied
#
varexp <- function(K=NULL, h2l=NULL, auc=NULL)
{
   if(is.null(K) || all(is.null(c(h2l, auc)))) {
      stop("Must specify K and at least one of h2l, auc")
   }

   if(any(c(K, h2l, auc) <= 0) || any(c(K, h2l, auc) > 1)) {
      stop("K, h2l, and auc must be in the range (0, 1]")
   }

   T <- qnorm(1 - K)
   z <- dnorm(T)
   i <- z / K
   v <- -i * K / (1 - K)
   
   aucmax <- as.numeric(NA)
   if(!is.null(h2l)) {
      # maximum achievable AUC
      mean1 <- i * h2l
      mean2 <- v * h2l
      var1 <- h2l * (1 - h2l * i * (i - T))
      var2 <- h2l * (1 - h2l * v * (v - T))
      
      d <- (mean1 - mean2) / sqrt(var1 + var2)
      aucmax <- pnorm(d)
   }
   
   h2lx <- rho2gg <- as.numeric(NA)
   if(!is.null(auc)) {
      # variance of genetic profile
      Q <- qnorm(auc, lower.tail=TRUE)
      h2lx <- 2 * Q^2 / ((v - i)^2 + Q^2 * i * (i - T) + v * (v - T))

      if(!is.null(h2l)) {
	 rho2gg <- h2lx / h2l
	 rho2gg <- pmin(rho2gg, 1)
      }
   }

   cbind(K=K, h2l=h2l, aucmax=aucmax, varexp=h2lx, genvarexp=rho2gg)
}

varexp.test <- function()
{
   K <- 0.0054
   h2l <- 0.86
   auc <- 0.96
   v <- varexp(K=K, h2l=h2l, auc=auc)
   aucmax <- 0.9930707
   h2lx <- 0.549721
   rho2gg <- 0.6392105

   if((aucmax - v[1])^2 > 1e-6)
      cat("bad aucmax:", aucmax, v[1], "\n")

   if((h2lx - v[2])^2 > 1e-6)
      cat("bad h2lx:", h2lx, v[2], "\n")
   
   if((rho2gg - v[3])^2 > 1e-6)
      cat("bad rho2gg:", rho2gg, v[3], "\n")

}

