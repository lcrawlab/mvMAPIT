library(doParallel)
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
library(CompQuadForm)

ggd.qqplot <- function(pvector, main = NULL, ...) {
  o <- -log10(sort(pvector, decreasing = F))
  e <- -log10(1:length(o) / length(o))
  plot(e, o,
    pch = 19, cex = 1, main = main, ...,
    xlab = expression(Expected ~ ~ -log[10](italic(p))),
    ylab = expression(Observed ~ ~ -log[10](italic(p))),
    xlim = c(0, max(e)), ylim = c(0, max(o))
  )
  lines(e, e, col = "red")
}

set.seed(11151990)


ind <- 3e2
nsnp <- 1e2
pve <- 0.6
rho <- 0.5
maf <- 0.05 + 0.45 * runif(nsnp)
X <- (runif(ind * nsnp) < maf) + (runif(ind * nsnp) < maf)
X <- matrix(as.double(X), ind, nsnp, byrow = TRUE)
Xmean <- apply(X, 2, mean)
Xsd <- apply(X, 2, sd)
X <- t((t(X) - Xmean) / Xsd)


ncausal1 <- 10
ncausal2 <- 10
ncausal3 <- 10



snp.ids <- 1:nsnp
s1 <- sample(snp.ids, ncausal1, replace = F)
s2 <- sample(snp.ids[-s1], ncausal2, replace = F)
s3 <- sample(snp.ids[-c(s1, s2)], ncausal3, replace = F)

Xcausal1 <- X[, s1]
Xcausal2 <- X[, s2]
Xcausal3 <- X[, s3]

W <- c()
for (i in 1:ncausal1) {
  W <- cbind(W, Xcausal1[, i] * Xcausal2)
}


Xmarginal <- cbind(Xcausal1, Xcausal2, Xcausal3)
beta <- rnorm(dim(Xmarginal)[2])
y_marginal <- Xmarginal %*% beta
beta <- beta * sqrt(pve * rho / var(y_marginal))
y_marginal <- Xmarginal %*% beta


alpha <- rnorm(dim(W)[2])
y_epi <- W %*% alpha
alpha <- alpha * sqrt(pve * (1 - rho) / var(y_epi))
y_epi <- W %*% alpha

y_err <- rnorm(ind)
y_err <- y_err * sqrt((1 - pve) / var(y_err))
y <- y_marginal + y_epi + y_err


colnames(X) <- paste("SNP", 1:ncol(X), sep = "")

SNPs <- colnames(X)[c(s1, s2)]


cores <- detectCores()


source("../R/MAPIT.R")
sourceCpp("../src/MAPIT.cpp")
# library(mvMAPIT)

ptm <- proc.time() # Start clock
# mapit = MAPIT(t(X),y,hybrid=FALSE,cores=cores)
hybrid <- MAPIT(t(X), y, hybrid = TRUE, cores = cores)
proc.time() - ptm # Stop clock

normal.pvals1 <- mapit$pvalues
names(normal.pvals1) <- colnames(X)
