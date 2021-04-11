#' Multivariate MArginal ePIstasis Test (mvMAPIT)
#' 
#' \code{MvMAPIT} will run a version of the MArginal ePIstasis Test (MAPIT) under the following model variations:
#' 
#' This function will run a multivariate version of the MArginal ePIstasis Test (mvMAPIT).
#' 
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
#' @param X is the p x n genotype matrix where p is the number of variants and n is the number of samples. Must be a matrix and not a data.frame.
#' @param y is the n x d matrix of quantitative or continuous traits.
#' @param W is the matrix qxn matrix of covariates. Must be a matrix and not a data.frame.
#' @param C is an nxn covariance matrix detailing environmental effects and population structure effects.
#' @param hybrid is a parameter detailing if the function should run the hybrid hypothesis testing procedure between the normal Z test and the Davies method. Default is TRUE.
#' @param threshold is a parameter detailing the value at which to recalibrate the Z test p values. If nothing is defined by the user, the default value will be 0.05 as recommended by the Crawford et al. (2017).
#' @param test is a parameter defining what hypothesis test should be implemented. Takes on values 'normal' or 'davies'. This parameter only matters when hybrid = FALSE. If test is not defined when hybrid = FALSE, the function will automatically use test = 'normal'.
#' @param cores is a parameter detailing the number of cores to parallelize over. It is important to note that this value only matters when the user has implemented OPENMP on their operating system. If OPENMP is not installed, then please leave cores = 1 and use the standard version of this code and software.
#' 
#' @return A list of P values and PVEs
#' @useDynLib mvMAPIT
#' @export
#' @import CompQuadForm
MvMAPIT <- function(X, y, W = NULL, C = NULL, hybrid = TRUE, threshold = 0.05, test = "normal", cores = 1) {

  if (cores > 1) {
    if (cores > detectCores()) {
      warning("The number of cores you're setting is larger than detected cores!")
      cores <- detectCores()
    }
  }

  if (hybrid == TRUE) {
    vc.mod <- MAPITCpp(X, y, W, C, NULL, "normal", cores = cores, NULL) # Normal Z-Test
    pvals <- vc.mod$pvalues
    names(pvals) <- rownames(X)
    pves <- vc.mod$PVE
    names(pves) <- rownames(X)

    ind <- which(pvals <= threshold) # Find the indices of the p-values that are below the threshold

    vc.mod <- MAPITCpp(X, y, W, C, ind, "davies", cores = cores, NULL)
    davies.pvals <- davies_exact(vc.mod, X)
    pvals[ind] <- davies.pvals[ind]
  } else if (test == "normal") {
    vc.mod <- MAPITCpp(X, y, W, C, NULL, "normal", cores = cores, NULL)
    pvals <- vc.mod$pvalues
    names(pvals) <- rownames(X)
    pves <- vc.mod$PVE
    names(pves) <- rownames(X)
  } else {
    ind <- 1:nrow(X)
    vc.mod <- MAPITCpp(X, y, W, C, ind, "davies", cores = cores, NULL)
    davies.pvals <- davies_exact(vc.mod, X)
    pvals <- davies.pvals
    pves <- vc.mod$PVE
  }
  return(list("pvalues" = pvals, "pves" = pves))
}

# Runs the Davies portion of the hypothesis testing
davies_exact <- function(vc.mod, X) {
  ### Apply Davies Exact Method ###
  vc.ts <- vc.mod$Est
  names(vc.ts) <- rownames(X)

  davies.pvals <- c()
  for (i in 1:length(vc.ts)) {
    lambda <- sort(vc.mod$Eigenvalues[, i], decreasing = T)

    Davies_Method <- davies(vc.mod$Est[i], lambda = lambda, acc = 1e-8)
    davies.pvals[i] <- 2 * min(Davies_Method$Qq, 1 - Davies_Method$Qq)
    names(davies.pvals)[i] <- names(vc.ts[i])
  }
  return(davies.pvals)
}
