#' Multivariate MArginal ePIstasis Test (mvMAPIT)
#'
#' \code{MvMAPIT} will run a version of the MArginal ePIstasis Test (MAPIT) under the following model variations:
#'
#' This function will run a multivariate version of the MArginal ePIstasis Test (mvMAPIT).
#'
#' (1) Standard Model: y = m+g+e where m ~ MVN(0,omega^2K), g ~ MVN(0,sigma^2G), e ~ MVN(0,tau^2M).
#' Recall from Crawford et al. (2017) that m is the combined additive effects from all other variants,
#' and effectively represents the additive effect of the k-th variant under the polygenic background
#' of all other variants; K is the genetic relatedness matrix computed using
#' genotypes from all variants other than the k-th; g is the summation of all pairwise interaction
#' effects between the k-th variant and all other variants; G represents a relatedness matrix
#' computed based on pairwise interaction terms between the k-th variant and all other variants. Here,
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
#' @param Y is the d x n matrix of d quantitative or continuous traits for n samples.
#' @param Z is the matrix q x n matrix of covariates. Must be a matrix and not a data.frame.
#' @param C is an n x n covariance matrix detailing environmental effects and population structure effects.
#' @param threshold is a parameter detailing the value at which to recalibrate the Z test p values. If nothing is defined by the user, the default value will be 0.05 as recommended by the Crawford et al. (2017).
#' @param accuracy is a parameter setting the davies function numerical approximation accuracy. This parameter is not needed for the normal test. Smaller p-values than the accuracy will be zero.
#' @param test is a parameter defining what hypothesis test should be implemented. Takes on values 'normal', 'davies', and 'hybrid'. The 'hybrid' test runs first the 'normal' test and then the 'davies' test on the significant variants.
#' @param cores is a parameter detailing the number of cores to parallelize over. It is important to note that this value only matters when the user has implemented OPENMP on their operating system. If OPENMP is not installed, then please leave cores = 1 and use the standard version of this code and software.
#' @param variantIndex is a vector containing indices of variants to be included in the computation.
#' @param phenotypeCovariance is a string parameter defining how to model the covariance between phenotypes of effects. Possible values: 'identity', 'covariance', 'homogeneous', 'combinatorial'.
#' @param logLevel is a string parameter defining the log level for the logging package.
#' @param logFile is a string parameter defining the name of the log file for the logging output.
#'
#' @return A list of P values and PVEs
#' @useDynLib mvMAPIT
#' @export
#' @import CompQuadForm
#' @import Rcpp
MvMAPIT <- function(X,
                    Y,
                    Z = NULL,
                    C = NULL,
                    threshold = 0.05,
                    accuracy = 1e-8,
                    test = c('normal', 'davies', 'hybrid'),
                    cores = 1,
                    variantIndex = NULL,
                    phenotypeCovariance = c('identity', 'covariance', 'homogeneous', 'combinatorial'),
                    logLevel = 'WARN',
                    logFile = NULL) {

  test <- match.arg(test)
  phenotypeCovariance <- match.arg(phenotypeCovariance)
  if (cores > 1) {
    if (cores > detectCores()) {
      warning("The number of cores you're setting is larger than detected cores!")
      cores <- detectCores()
    }
  }

  logging::logReset()
  logging::basicConfig(level = logLevel)
  log <- logging::getLogger('MvMAPIT')
  if(!is.null(logFile)) {
    filePath <- file.path(getwd(),logFile)
    log$debug('Logging to file: %s', filePath)
    log$addHandler(logging::writeToFile, file=filePath)
  }

  if (is.vector(Y)) {
    Y <- t(Y)
  }

  log$debug('Running in %s test mode.', test)
  log$debug('Phenotype covariance: %s', phenotypeCovariance)
  log$debug('Genotype matrix: %d x %d', nrow(X), ncol(X))
  log$debug('Phenotype matrix: %d x %d', nrow(Y), ncol(Y))
  log$debug('Genotype matrix determinant: %f', det((X) %*% t(X)))
  zero_var <- which(apply(X, 1, var) == 0)
  log$debug('Number of zero variance variants: %d', length(zero_var))
  X <- remove_zero_variance(X) # operates on rows
  log$debug('Genotype matrix after removing zero variance variants: %d x %d', nrow(X), ncol(X))

  if (test == 'hybrid') {
    vc.mod <- MAPITCpp(X, Y, Z, C, variantIndex, "normal", cores = cores, NULL, phenotypeCovariance) # Normal Z-Test
    pvals <- vc.mod$pvalues
    #row.names(pvals) <- rownames(X)
    pves <- vc.mod$PVE
    #row.names(pves) <- rownames(X)
    timings <- vc.mod$timings
    ind <- which(pvals <= threshold) # Find the indices of the p-values that are below the threshold
    if (phenotypeCovariance == 'combinatorial') {
      any_significance <- apply(pvals, 1, function(r) any(r <= threshold))
      ind_temp <- ind
      ind <- which(any_significance == TRUE)
    }
    log$info('%d p-values are significant with alpha = %f', length(ind), threshold)

    log$info('Running davies method on selected SNPs.')
    vc.mod <- MAPITCpp(X, Y, Z, C, ind, "davies", cores = cores, NULL, phenotypeCovariance)
    davies.pvals <- mvmapit_pvalues(vc.mod, X, accuracy)
    if (phenotypeCovariance == 'combinatorial') {
      pvals[ind_temp] <- davies.pvals[ind_temp]
    } else {
      pvals[ind] <- davies.pvals[ind]
    }
  } else if (test == "normal") {
    vc.mod <- MAPITCpp(X, Y, Z, C, variantIndex, "normal", cores = cores, NULL, phenotypeCovariance)
    pvals <- vc.mod$pvalues
    pves <- vc.mod$PVE
    timings <- vc.mod$timings
  } else {
    ind <- ifelse(variantIndex, variantIndex, 1:nrow(X))
    vc.mod <- MAPITCpp(X, Y, Z, C, ind, "davies", cores = cores, NULL, phenotypeCovariance)
    pvals <- mvmapit_pvalues(vc.mod, X, accuracy)
    pves <- vc.mod$PVE
    timings <- vc.mod$timings
  }
  timings_mean <- apply(as.matrix(timings[rowSums(timings) != 0, ]), 2, mean)
  log$info('Calculated mean time of execution. Return list.')
  row.names(pvals) <- rownames(X)
  row.names(pves) <- rownames(X)
  if (length(rownames(Y)) > 0) {
    column_names <- mapit_struct_names(Y, phenotypeCovariance)
  } else if (nrow(Y) > 1 && (phenotypeCovariance == 'combinatorial')) {
    row.names(Y) <- sprintf("P%s", 1:nrow(Y))
    column_names <- mapit_struct_names(Y, phenotypeCovariance)
  } else {
    column_names <- NULL
  }
  colnames(pvals) <- column_names
  colnames(pves) <- column_names
  if (!is.null(variantIndex)) {
    log$debug('Set pve to NA if not in varianIndex.')
    pves[!(c(1:nrow(pves)) %in% variantIndex)] <- NA
  }
  if (nrow(Y) > 1 && (phenotypeCovariance == 'combinatorial')) {
      pves <- set_covariance_pves(Y, pves)
  }
  return(list("pvalues" = pvals, "pves" = pves, "timings" = timings_mean))
}

remove_zero_variance <- function(X) {
  return(X[which(apply(X, 1, var) != 0),])
}

# This naming sequence has to match the creation of the q-matrix in the C++ routine of mvMAPIT
mapit_struct_names <- function (Y, phenotypeCovariance) {
  if (length(phenotypeCovariance) > 0 && !(phenotypeCovariance == 'combinatorial')) {
    return(c('kronecker'))
  }
  phenotype_names <- rownames(Y)
  phenotype_combinations <- c()
  for (i in seq_len(nrow(Y))) {
    for (j in seq_len(nrow(Y))) {
      if (j <= i) {
        phenotype_combinations <- c(phenotype_combinations,
                                    paste0(phenotype_names[i], "*", phenotype_names[j]))
      }
    }
  }
  return(phenotype_combinations)
}

set_covariance_pves  <- function(Y, pves) {
  counter <- 0
  for (i in seq_len(nrow(Y))) {
    for (j in seq_len(nrow(Y))) {
      if (j <= i) {
        counter <- counter + 1
        if (j < i) pves[, counter] <- NA
      }
    }
  }
  return(pves)
}
