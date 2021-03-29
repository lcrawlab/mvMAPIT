#' simulate
#' 
#' Simulate data. Deprecated.
#' 
#' This function will run a version of the MArginal ePIstasis Test (MAPIT) under the following model variations:
#' (1) Standard Model: y = m+g+e where m ~ MVN(0,omega^2K), g ~ MVN(0,sigma^2G), e ~ MVN(0,tau^2M).
#' Recall from Crawford et al. (2017) that m is the combined additive effects from all other variants,
#' and effectively represents the additive effect of the kth variant under the polygenic background
#' of all other variants; K is the genetic relatedness matrix computed using
#' genotypes from all variants other than the kth; g is the summation of all pairwise interaction
#' effects between the kth variant and all other variants; G represents a relatedness matrix
#' computed based on pairwise interaction terms between the kth variant and all other variants. Here,
#' we also denote D = diag(x_k) to be an n × n diagonal matrix with the genotype vector x_k as its
#' diagonal elements. It is important to note that both K and G change with every new marker k that is
#' considered. Lastly; M is a variant specific projection matrix onto both the null space of the intercept
#' and the corresponding genotypic vector x_k.
#' 
#' (2) Standard + Covariate Model: y = Wa+m+g+e where W is a matrix of covariates with effect sizes a.
#' 
#' (3) Standard + Common Environment Model: y = m+g+c+e where c ~ MVN(0,eta^2C) controls for extra
#' environmental effects and population structure with covariance matrix C.
#' 
#' (4) Standard + Covariate + Common Environment Model: y = Wa+m+g+c+e
#' 
#' This function will consider the following three hypothesis testing strategies which are featured in Crawford et al. (2017):
#' (1) The Normal or Z test
#' (2) Davies Method
#' (3) Hybrid Method (Z test + Davies Method)
#' 
#' @param nsnps Number of causal SNPs
#' @param npepi Number of pleiotropic epistatic SNPs. npepi+ntepi must be less than nsnps.
#' @param ntepi Number of trait specific epistatic SNPs. npepi+ntepi must be less than nsnps.
#' @param H2 Broad-sense heritability. Default 0.6.
#' @param d Number of phenotypes. Default 2.
#' @param rho Proportion of heritability explained by additivity. Default 0.8.
#' @param seed Random seed for simulation. Default 67132.
#' @param data Data frame from which to draw SNPs. Default random_genotype_matrix.
#' @param outfile is the filename for the output file.
#' @return A list of P values and PVEs
#' @useDynLib mvMAPIT
#' @export
#' @import CompQuadForm
#' @import doParallel
#' @import Rcpp
#' @import RcppArmadillo
#' @import Matrix
#' @import mvtnorm
#' @import PHENIX
simulate <- function(nsnps,
                     npepi,
                     ntepi,
                     H2 = 0.6,
                     d = 2,
                     rho = 0.8,
                     seed = 67132,
                     data = "random_genotype_matrix",
                     outfile = NULL) {
  set.seed(seed)
  load.Rdata(data)

  # subset data
  random_genotype_matrix <- random_genotype_matrix[1:nsnps, 1:nsnps]
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

  npadd <- (nsnps - (npepi + ntepi)) / 2 # pleiotropic additive
  ntadd <- (nsnps - (npepi + ntepi)) / 2 # trait specific additive

  stopifnot(npadd > 0, ntadd > 0) # check that npadd and ntadd are positive

  ### Set the Correlations ###
  add.cor <- 0
  epi.cor <- 0
  env.cor <- 0

  ### Simulate the Data ###
  Yn1 <- c()
  Yn2 <- c()
  S1 <- c()


  s <- 1:nsnp
  s1 <- sample(s, npepi, replace = F)
  s3 <- sample(s[s %in% s1 == FALSE], npadd, replace = F)

  for (j in 1:d) {
    # Select Causal SNPs
    s2 <- sample(s[s %in% c(s1, s3) == FALSE], ntepi / 2, replace = F)
    s4 <- sample(s[s %in% c(s1, s2, s3) == FALSE], ntadd / 2, replace = F)

    # Create Causal Epistatic Matrix
    Xcausal1 <- X[, s1]
    Xcausal2 <- X[, s2]
    Xcausal3 <- X[, s3]
    Xcausal4 <- X[, s4]
    Xepi <- c()
    for (i in 1:npepi) {
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
  # Ymean=apply(Y, 2, mean); Ysd=apply(Y, 2, sd); Y=t((t(Y)-Ymean)/Ysd)

  return(s1, s2, s3, s4, X, Y)
}

### Scenario Combinations {p1/p2}: (I) 10/20; (II) 10/100; (III) 20/20; (IV) 20/100 ###



# # run MAPIT
# library('Rcpp')
# source('R/MAPIT.R')
# sourceCpp('../src/MAPIT.cpp')
# cores = detectCores()
#
# # if hybrid=FALSE, should return same results as regular MAPIT, yes?
#
# ptm <- proc.time() #Start clock
# mapit = MvMAPIT(t(X), Y[,1], Y[,2], hybrid=FALSE, cores=cores)
# proc.time() - ptm #Stop clock
#
# normal.pvals1 = mapit$pvalues
# names(normal.pvals1) = colnames(X)
# save(mapit, file="mapit.Rdata")
#
# normal.pvals1[s1] # pvalues for pleiotropic epistatic
#
# ptm <- proc.time() #Start clock
# hybrid = MvMAPIT(t(X), Y[,1], Y[,2], hybrid=TRUE, cores=cores)
# proc.time() - ptm #Stop clock
#
# hybrid.pvals = hybrid$pvalues
# names(hybrid.pvals) = colnames(X)
# save(hybrid, file="hybrid.Rdata")
#####################################################################################
######################################################################################
######################################################################################
