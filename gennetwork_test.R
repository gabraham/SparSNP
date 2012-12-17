
set.seed(2879)

Y <- round(matrix(rnorm(20 * 5), 20, 5), 2)

cortype <- 1
corthresh <- 0

R <- cor(Y)

R[abs(R) < corthresh] <- 0

K <- ncol(Y)
nV <- K
UR <- R
UR[lower.tri(UR)] <- 0
diag(UR) <- 0

W <- if(cortype == 0) {
   (abs(UR) > corthresh) + 0
} else if(cortype == 1) {
   abs(UR)
} else if(cortype == 2) {
   UR^2
} else {
   stop("unknown cortype:", cortype)
}

nzUR <- which(UR != 0)

E <- which(UR != 0, arr.ind=TRUE)
Ecoef <- W[nzUR]
Esign <- sign(R[nzUR])

#nE <- nrow(E)
nE <- K * (K - 1) / 2

C_I <- c(1:nE, 1:nE)
C_J <- as.numeric(E)
C_S <- cbind(Ecoef, -Ecoef * Esign)
C <- matrix(0, nE, nV)
C[cbind(C_I, C_J)] <- C_S

system("./gennetwork_test")
C2 <- matrix(scan("C.txt", sep=","), nE, nV, byrow=TRUE)

mean((C2 - C)^2)


