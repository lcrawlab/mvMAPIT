### Clear Console ###
cat("\014")

### Clear Environment ###
rm(list = ls(all = TRUE))

### Load in the R libraries ###
library(doParallel)
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
library(CompQuadForm)
library(Matrix)
library(mvtnorm)

### Load in the Phenix Package ###
library(PHENIX) # (From the Marchini Group)

######################################################################################
######################################################################################
######################################################################################

### Set seed ###
set.seed(67132)

### Load in the Data ###
# load("random_genotype_matrix.rda")
# subset data to run faster
# random_genotype_matrix = random_genotype_matrix[1:1000, 1:1000] #ok, because here we just do 500
random_genotype_matrix <- random_genotype_matrix[1:500, 1:500]
# if we comment out the line it works
# random_genotype_matrix = random_genotype_matrix

unique.snps <- apply(random_genotype_matrix, 1, function(x) {
  length(unique(x))
})
random_genotype_matrix <- random_genotype_matrix[unique.snps > 1, ]
maf <- apply(random_genotype_matrix, 1, mean) / 2
X <- t(random_genotype_matrix[maf > 0.05, ])

Xmean <- apply(X, 2, mean)
Xsd <- apply(X, 2, sd)
X <- t((t(X) - Xmean) / Xsd)
ind <- dim(X)[1]
nsnp <- dim(X)[2]

H2 <- 0.6 # H2 = Broad-sense heritability
d <- 2 # d = Number of phenotypes
rho <- 0.8 # rho = Proportion of heritability explained by additivity
ncausal1 <- 10 # Set 1 of causal pleiotropic set of epistatic SNPs
ncausal2 <- 20 # Set 2 of causal trait specific set of epistatic SNPs

### Scenario Combinations {p1/p2}: (I) 10/20; (II) 10/100; (III) 20/20; (IV) 20/100 ###

### With these parameters, the combined interaction effect size in Scenarios I and II should be cut in half for Scenarios III and IV ###

# (I) ncausal1 = 10, ncausal2 = 20
# (II) ncausal1 = 10, ncausal2 = 100
# (III) ncausal1 = 20, ncausal2 = 20
# (IV) ncausal1 = 20, ncausal2 = 100

# ncausal3 = (1e3-(ncausal1+ncausal2))/2 #Set 3 of causal pleiotropic set of additive SNPs
# ncausal4 = (1e3-(ncausal1+ncausal2))/2 #Set 4 of trait specific set of additive SNPs
ncausal3 <- (500 - (ncausal1 + ncausal2)) / 2
ncausal4 <- (500 - (ncausal1 + ncausal2)) / 2

### Set the Correlations ###
add.cor <- 0
epi.cor <- 0
env.cor <- 0

### Simulate the Data ###
Yn1 <- c()
Yn2 <- c()
S1 <- c()


s <- 1:nsnp
s1 <- sample(s, ncausal1, replace = F)
s3 <- sample(s[s %in% s1 == FALSE], ncausal3, replace = F) # this is where the error is
# ok, because we are sampling from a set of less than 500

for (j in 1:d) {
  ## Select Causal SNPs
  s2 <- sample(s[s %in% c(s1, s3) == FALSE], ncausal2 / 2, replace = F)
  s4 <- sample(s[s %in% c(s1, s2, s3) == FALSE], ncausal4 / 2,
    replace =
      F
  )

  # Create Causal Epistatic Matrix
  Xcausal1 <- X[, s1]
  Xcausal2 <- X[, s2]
  Xcausal3 <- X[, s3]
  Xcausal4 <- X[, s4]
  Xepi <- c()
  for (i in 1:ncausal1) {
    Xepi <- cbind(Xepi, Xcausal1[, i] * Xcausal2)
  }
  dim(Xepi)

  # Marginal Effects Only
  Xmarginal <- X[, c(s1, s2, s3, s4)]
  beta <- rnorm(dim(Xmarginal)[2])
  y_marginal <- c(Xmarginal %*% beta)

  # Pairwise Epistatic Effects
  beta <- rnorm(dim(Xepi)[2])
  y_epi <- c(Xepi %*% beta)

  Yn1 <- cbind(Yn1, y_marginal)
  Yn2 <- cbind(Yn2, y_epi)
  S1 <- cbind(S1, s1)
}

save(s1, file = "s1.Rdata")
save(s2, file = "s2.Rdata")
save(s3, file = "s3.Rdata")
save(s4, file = "s4.Rdata")

### Set Additive Genetic Correlation Across Traits ###
B <- diag(d)
B <- diag(sqrt(H2 * rho), d) %*% cov2cor(B) %*% diag(sqrt(H2 * rho), d)

corr <- add.cor
C <- matrix(corr, ncol = d, nrow = d)
diag(C) <- 0
C <- C + diag(d) # Environmental Correlation Matrix
Yn1 <- Yn1 %*% mat.sqrt(C)
Yn1 <- Yn1 %*% mat.sqrt(B / apply(Yn1, 2, var))

# cov(Yn1); cor(Yn1)

### Set Epistatic Genetic Correlation Across Traits ###
B <- diag(d)
B <- diag(sqrt(H2 * (1 - rho)), d) %*% cov2cor(B) %*% diag(sqrt(H2 * (1 - rho)), d)

corr <- epi.cor
C <- matrix(epi.cor, ncol = d, nrow = d)
diag(C) <- 0
C <- C + diag(d) # Environmental Correlation Matrix
Yn2 <- Yn2 %*% mat.sqrt(C)
Yn2 <- Yn2 %*% mat.sqrt(B / apply(Yn2, 2, var))

# cov(Yn2); cor(Yn2)

### Do the Environmental Effects ###
V <- diag(d) # Environmental Covariance Matrix
V <- diag(sqrt(1 - H2), d) %*% cov2cor(V) %*% diag(sqrt(1 - H2), d)

corr <- env.cor
C <- matrix(corr, ncol = d, nrow = d)
diag(C) <- 0
C <- C + diag(d) # Environmental Correlation Matrix
E <- matrix(rnorm(ind * d), ind, d) %*% mat.sqrt(C)
E <- E %*% mat.sqrt(V / apply(E, 2, var))

# cov(E); cor(E)

### Compute the Phenotypes ###
Y <- Yn1 + Yn2 + E

### Check dimensions and add SNP names ###
dim(X)
dim(Y)
colnames(X) <- paste("SNP", 1:ncol(X), sep = "")
SNPs <- colnames(X)[s1]
save(X, file = "X.Rdata")
save(Y, file = "Y.Rdata")
# Ymean=apply(Y, 2, mean); Ysd=apply(Y, 2, sd); Y=t((t(Y)-Ymean)/Ysd)

# run MAPIT
source("../R/MAPIT.R")
sourceCpp("../src/MAPIT.cpp")
cores <- detectCores()

# if hybrid=FALSE, should return same results as regular MAPIT, yes?

ptm <- proc.time() # Start clock
mapit <- MvMAPIT(t(X), Y[, 1], Y[, 2], hybrid = FALSE, cores = cores)
proc.time() - ptm # Stop clock

normal.pvals1 <- mapit$pvalues
names(normal.pvals1) <- colnames(X)
save(mapit, file = "mapit.Rdata")

normal.pvals1[s1] # pvalues for pleiotropic epistatic

ptm <- proc.time() # Start clock
hybrid <- MvMAPIT(t(X), Y[, 1], Y[, 2], hybrid = TRUE, cores = cores)
proc.time() - ptm # Stop clock

hybrid.pvals <- hybrid$pvalues
names(hybrid.pvals) <- colnames(X)
save(hybrid, file = "hybrid.Rdata")
#####################################################################################
######################################################################################
######################################################################################
