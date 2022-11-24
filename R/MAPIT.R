#' Multivariate MArginal ePIstasis Test (mvMAPIT)
#'
#' This function runs a multivariate version of the MArginal ePIstasis
#' Test (mvMAPIT) under the following model variations:
#'
#' (1) Standard Model: y = m+g+e
#' where m ~ MVN(0,omega^2K), g ~ MVN(0,sigma^2G), e ~ MVN(0,tau^2M).
#' Recall from Crawford et al. (2017) that m is the combined additive effects
#' from all other variants, represents the additive effect of the k-th variant
#' under the polygenic background of all other variants; K is the genetic
#' relatedness matrix computed using genotypes from all variants other than the
#' k-th; g is the summation of all pairwise interaction effects between the
#' k-th variant and all other variants; G represents a relatedness matrix
#' computed based on pairwise interaction terms between the k-th variant and all
#' other variants. Here, we also denote D = diag(x_k) to be an n × n diagonal
#' matrix with the genotype vector x_k as its diagonal elements. It is important
#' to note that both K and G change with every new marker k that is considered.
#' Lastly; M is a variant specific projection matrix onto both the null space of
#' the intercept and the corresponding genotypic vector x_k.
#'
#' (2) Standard + Covariate Model: y = Wa+m+g+e
#' where W is a matrix of covariates with effect sizes a.
#'
#' (3) Standard + Common Environment Model: y = m+g+c+e i
#' where c ~ MVN(0,eta^2C) controls for extra environmental effects and
#' population structure with covariance matrix C.
#'
#' (4) Standard + Covariate + Common Environment Model: y = Wa+m+g+c+e
#'
#' This function will consider the following three hypothesis testing strategies
#' which are featured in Crawford et al. (2017):
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
#' @param test is a parameter defining what hypothesis test should be run. Takes on values 'normal', 'davies', and 'hybrid'. The 'hybrid' test runs first the 'normal' test and then the 'davies' test on the significant variants.
#' @param cores is a parameter detailing the number of cores to parallelize over. It is important to note that this value only matters when the user has installed OPENMP on their operating system.
#' @param variantIndex is a vector containing indices of variants to be included in the computation.
#' @param logLevel is a string parameter defining the log level for the logging package.
#' @param logFile is a string parameter defining the name of the log file for the logging output. Default is stdout.
#'
#' @return A list of P values and PVEs
#' @examples
#' set.seed(837)
#' p <- 200
#' n <- 100
#' d <- 2
#' X <- matrix(
#'     runif(p * n),
#'     ncol = p
#' )
#' Y <- matrix(
#'     runif(d * n),
#'     ncol = d
#' )
#' mapit <- mvmapit(
#'     t(X),
#'     t(Y),
#'     test = "normal", cores = 1, logLevel = "INFO"
#' )
#' @useDynLib mvMAPIT
#' @name mvmapit
#' @export
#' @import CompQuadForm
#' @import Rcpp
mvmapit <- function(
    X, Y, Z = NULL, C = NULL, threshold = 0.05, accuracy = 1e-08, test = c("normal", "davies", "hybrid"),
    cores = 1, variantIndex = NULL, logLevel = "WARN", logFile = NULL
) {

    test <- match.arg(test)
    if (cores > 1) {
        if (cores > detectCores()) {
            warning("The number of cores you're setting is larger than detected cores!")
            cores <- detectCores()
        }
    }

    logging::logReset()
    logging::basicConfig(level = logLevel)
    log <- logging::getLogger("mvmapit")
    if (!is.null(logFile)) {
        filePath <- file.path(getwd(), logFile)
        log$debug("Logging to file: %s", filePath)
        log$addHandler(logging::writeToFile, file = filePath)
    }

    if (is.vector(Y)) {
        Y <- t(Y)
    }
    row.names(Y) <- make.unique(as.character(row.names(Y)))

    log$debug("Running in %s test mode.", test)
    log$debug(
        "Genotype matrix: %d x %d", nrow(X),
        ncol(X)
    )
    log$debug(
        "Phenotype matrix: %d x %d", nrow(Y),
        ncol(Y)
    )
    zero_var <- which(apply(X, 1, var) == 0)
    log$debug("Number of zero variance variants: %d", length(zero_var))
    X <- remove_zero_variance(X)  # operates on rows
    log$debug(
        "Genotype matrix after removing zero variance variants: %d x %d", nrow(X),
        ncol(X)
    )

    log$debug("Scale X matrix appropriately.")
    Xsd <- apply(X, 1, sd)
    Xmean <- apply(X, 1, mean)
    X <- (X - Xmean) / Xsd

    variance_components_ind <- get_variance_components_ind(Y)
    if (test == "hybrid") {
        vc.mod <- MAPITCpp(X, Y, Z, C, variantIndex, "normal", cores = cores, NULL)  # Normal Z-Test
        pvals <- vc.mod$pvalues
        # row.names(pvals) <- rownames(X)
        pves <- vc.mod$PVE
        # row.names(pves) <- rownames(X)
        timings <- vc.mod$timings
        ind_matrix <- which(pvals[, variance_components_ind] <= threshold)  # Find the variance component indices of the p-values that are below the threshold
        log$info(
            "%d p-values of variance components are significant with alpha = %f",
            length(ind_matrix),
            threshold
        )
        if (nrow(Y)) {
            any_significance <- apply(pvals, 1, function(r) any(r <= threshold))
        } else {
            any_significance <- apply(pvals[, variance_components_ind], 1, function(r) any(r <= threshold))
        }
        ind <- which(any_significance == TRUE)
        log$info(
            "%d positions are significant with alpha = %f", length(ind),
            threshold
        )

        log$info("Running davies method on selected SNPs.")
        vc.mod <- MAPITCpp(X, Y, Z, C, ind, "davies", cores = cores, NULL)
        davies.pvals <- mvmapit_pvalues(vc.mod, X, accuracy)
        pvals[, variance_components_ind][ind_matrix] <- davies.pvals[, variance_components_ind][ind_matrix]
    } else if (test == "normal") {
        vc.mod <- MAPITCpp(X, Y, Z, C, variantIndex, "normal", cores = cores, NULL)
        pvals <- vc.mod$pvalues
        pves <- vc.mod$PVE
        timings <- vc.mod$timings
    } else {
        ind <- ifelse(variantIndex, variantIndex, 1:nrow(X))
        vc.mod <- MAPITCpp(X, Y, Z, C, ind, "davies", cores = cores, NULL)
        pvals <- mvmapit_pvalues(vc.mod, X, accuracy)
        pvals <- set_covariance_components(variance_components_ind, pvals)
        pves <- vc.mod$PVE
        timings <- vc.mod$timings
    }
    timings_mean <- apply(as.matrix(timings[rowSums(timings) != 0, ]), 2, mean)
    log$debug("Calculated mean time of execution. Return list.")
    row.names(pvals) <- rownames(X)
    row.names(pves) <- rownames(X)
    column_names <- mapit_struct_names(Y)
    colnames(pvals) <- column_names
    colnames(pves) <- column_names
    if (!is.null(variantIndex)) {
        log$debug("Set pve to NA if not in varianIndex.")
        pves[!(c(1:nrow(pves)) %in%
            variantIndex)] <- NA
        pvals[!(c(1:nrow(pvals)) %in%
            variantIndex)] <- NA
    }
    pvals_df <- as.data.frame(pvals)
    pvals <- pvals_df %>%
        mutate(id = row.names(pvals_df)) %>%
        tidyr::pivot_longer(cols = !id,
                     names_to = "trait", values_to = "p")
    pves_df <- as.data.frame(pves)
    pves <- pves_df %>%
        mutate(id = row.names(pves_df)) %>%
        tidyr::pivot_longer(cols = !id,
                     names_to = "trait", values_to = "PVE")
    duration_ms <- timings_mean
    process <- c("cov", "projections", "vectorize", "q", "S", "vc")
    timings_mean <- data.frame(process, duration_ms)
    return(list(pvalues = pvals, pves = pves, duration = timings_mean))
}

#' Remove variants that don't vary accross the genotype data.
#'
#' This function takes in the genotype matrix and reomces the varaints with zero
#' variance.
#'
#' @param X Genotype matrix.
#' @return Genotype matrix X without zero variance variants.
#' @noRd
remove_zero_variance <- function(X) {
    return(X[which(apply(X, 1, var) != 0), ])
}

#' Create names for the columns of the mvMAPIT return object.
#'
#' This function takes in the traits matrix and creates names for the p-values
#' and PVEs for all variance and covariance components.
#' This naming sequence has to match the creation of the q-matrix in the C++
#' routine of mvMAPIT.
#'
#' @param Y Trait matrix.
#' @return Vector of strings for all trait combinations.
#' @noRd
mapit_struct_names <- function(Y) {
    phenotype_names <- rownames(Y)
    if (length(phenotype_names) == 0) {
        phenotype_names <- sprintf("P%s", 1:nrow(Y))
    }
    if (length(phenotype_names) == 1) {
        return(phenotype_names)
    }
    phenotype_combinations <- c()
    for (i in seq_len(nrow(Y))) {
        for (j in seq_len(nrow(Y))) {
            if (j <= i) {
                phenotype_combinations <- c(phenotype_combinations,
                                            sprintf(
                                                "%s*%s",
                                                phenotype_names[i],
                                                phenotype_names[j]
                                            )
                )
            }
        }
    }
    return(phenotype_combinations)
}

#' Set covariance components as "NA".
#'
#' This function takes in the column indices for the covariance components and
#' a matrix with data for both variance and covariance components and sets the
#' covariance components to NA. This is needed for the davies method version of
#' the p-values computation since there is no Davies method implementation for
#' p-values currently.
#'
#' @param variance_components_ind Column index vector for variance components.
#' @param X Data with both variance and covariance components.
#' @return Vector of strings for all trait combinations.
#' @noRd
set_covariance_components <- function(variance_components_ind, X) {
    X[, !variance_components_ind] <- NA
    return(X)
}

#' Compute index vector for columns containing data for variance components.
#'
#' This function takes in the trait matrix and computes the column indices that
#' contain the variance component portion in the return data from the C++
#' method.
#'
#' @param Y Trait matrix.
#' @return Vector of indices for all variance components.
#' @noRd
get_variance_components_ind <- function(Y) {
    ind <- c()
    counter <- 0
    for (i in seq_len(nrow(Y))) {
        for (j in seq_len(nrow(Y))) {
            if (j <= i) {
                counter <- counter + 1
                if (j == i) {
                  ind[counter] <- TRUE
                } else {
                  ind[counter] <- FALSE
                }
            }
        }
    }
    return(ind)
}

