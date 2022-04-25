#' Simulate phenotye data
#'
#' \code{simulate_phenotype} simulates phenotype data from a genotype matrix.
#'
#' This function takes a genotype matrix and simulates phenotype data under the following model:
#'
#' beta_i ~ MN(0, V_i, I), i in \{ additive, epistatic, residual\}
#'
#' The effect sizes follow a matrix normal distribution with no correlation between the samples but covariance between the effects for different phenotypes
#'
#'
#'
#'
#' @param genotype_matrix Genotype matrix with samples as rows, and SNPs as columns.
#' @param n_causal Number of SNPs that are causal.
#' @param n_trait_specific Number of causal SNPs with single trait epistatic effects.
#' @param n_pleiotropic Number of SNPs with pleiotropic effects.
#' @param group_ratio_trait Ratio of sizes of trait specific groups that interact, e.g. a ratio 1:3 would be value 3.
#' @param group_ratio_pleiotropic Ratio of sizes of pleiotropic groups that interact, e.g. a ratio 1:3 would be value 3.
#' @param H2 Broad-sense heritability. Can be vector.
#' @param d Number of phenotypes.
#' @param rho Proportion of heritability explained by additivity.
#' @param marginal_correlation Correlation between the additive effects of the phenotype.
#' @param epistatic_correlation Correlation between the epistatic effects of the phenotype.
#' @param seed Random seed for simulation.
#' @param logLevel is a string parameter defining the log level for the logging package.
#' @param logFile is a string parameter defining the name of the log file for the logging output.
#' @param maf_threshold is a float parameter defining the threshold for the minor allele frequency not included in causal SNPs.
#' @return A list object containing the phenotype data, the genotype data, as well as the causal SNPs and summary statistics.
#' @useDynLib mvMAPIT
#' @export
#' @import checkmate
#' @import dplyr
#' @import foreach
#' @import mvtnorm
#' @import parallel
simulate_phenotypes <- function(
    genotype_matrix, n_causal = 1000, n_trait_specific = 10, n_pleiotropic = 10,
    H2 = 0.6, d = 2, rho = 0.8, marginal_correlation = 0.3, epistatic_correlation = 0.3,
    group_ratio_trait = 1, group_ratio_pleiotropic = 1, maf_threshold = 0.01, seed = 67132,
    logLevel = "INFO", logFile = NULL
) {

    heritability <- rep(H2, d)[1:d]
    set.seed(seed)
    coll <- makeAssertCollection()
    assertInt(n_causal, lower = 0, add = coll)
    assertInt(n_trait_specific, lower = 0, add = coll)
    assertInt(n_pleiotropic, lower = 0, add = coll)
    assertDouble(group_ratio_trait, lower = 1, add = coll)
    assertDouble(group_ratio_pleiotropic, lower = 1, add = coll)
    assertDouble(heritability, lower = 0, upper = 1, add = coll)
    assertDouble(rho, lower = 0, upper = 1, add = coll)
    assertDouble(marginal_correlation, lower = -1, upper = 1, add = coll)
    assertDouble(epistatic_correlation, lower = -1, upper = 1, add = coll)
    assertDouble(maf_threshold, lower = 0, upper = 1, add = coll)
    assertInt(d, lower = 1, add = coll)
    assertInt(seed, lower = 1, add = coll)
    assertMatrix(genotype_matrix, all.missing = FALSE, add = coll)
    reportAssertions(coll)

    logging::logReset()
    logging::basicConfig(level = logLevel)
    log <- logging::getLogger("simulate_phenotypes")
    if (!is.null(logFile)) {
        filePath <- file.path(getwd(), logFile)
        log$debug("Logging to file: %s", filePath)
        log$addHandler(logging::writeToFile, file = filePath)
    }

    snp.ids <- 1:ncol(genotype_matrix)
    maf <- colMeans(genotype_matrix) / 2
    X <- scale(genotype_matrix)
    maf_compliant <- (maf > maf_threshold) & (maf < 1 - maf_threshold)
    # scale produces NaN when the columns have zero variance
    snp.ids.filtered <- snp.ids[complete.cases(t(X)) & maf_compliant]

    n_samples <- nrow(X)  # number of genotype samples
    n_snp <- length(snp.ids.filtered)  # number of SNPs passing quality control
    log$debug("Scaled genotype matrix: %d x %d", n_samples, n_snp)
    log$debug(
        "Disregard %d variants due to zero variance or small minor allele frequency.",
        ncol(genotype_matrix) - length(snp.ids.filtered)
    )
    log$debug("Minor allele frequency threshold %f.", maf_threshold)

    # divide groups into ratios
    n_group1_trait = ceiling(n_trait_specific / (1 + group_ratio_trait))
    n_group1_pleiotropic = ceiling(n_pleiotropic / (1 + group_ratio_pleiotropic))

    coll <- makeAssertCollection()
    assertInt(n_causal, lower = 0, upper = n_snp, add = coll)
    assertInt(n_pleiotropic + n_trait_specific, lower = 0, upper = n_causal, add = coll)
    assertInt(n_group1_trait, lower = 0, upper = n_trait_specific, add = coll)
    assertInt(n_group1_pleiotropic, lower = 0, upper = n_pleiotropic, add = coll)
    reportAssertions(coll)

    # factor vectors for splitting the groups
    f_trait <- get_factors(n_group1_trait, n_trait_specific)
    f_pleiotropic <- get_factors(n_group1_pleiotropic, n_pleiotropic)


    log$debug("Number of causal SNPs: %d", n_causal)
    log$debug("Number of trait specific SNPs: %d", n_trait_specific)
    log$debug("Number of pleiotropic SNPs: %d", n_pleiotropic)
    log$debug("NA in raw genotype matrix: %d", sum(is.na(genotype_matrix)))
    log$debug("NA in scaled genotype matrix: %d", sum(is.na(X)))

    Y <- c()
    causal_snps <- list()

    pleiotropic_set <- sample(snp.ids.filtered, n_pleiotropic, replace = F)  # declare peleiotropic SNPs before since they have to be present in every phenotype
    pleio_split <- split(pleiotropic_set, f = f_pleiotropic)
    X_pleio_group1 <- X[, pleio_split$group1]
    X_pleio_group2 <- X[, pleio_split$group2]
    X_epi_pleio <- foreach(i = seq_len(n_group1_pleiotropic), .combine = cbind) %do% {
            # this step fails if there are too little pleiotropic SNPs; i.e. <=
            # 3?
            X_pleio_group1[, i] * X_pleio_group2
    }
    if (!is.matrix(X_epi_pleio)) {
        X_epi_pleio <- matrix(
            NA, nrow = nrow(X),
            ncol = 0
        )
    }
    log$debug(
        "Dimensions of pleiotropic interaction matrix: %s x %s", nrow(X_epi_pleio),
        ncol(X_epi_pleio)
    )

    log$debug("Draw effects from multivariate normal with desired correlation.")

    C_marginal <- matrix(marginal_correlation, ncol = d, nrow = d)
    diag(C_marginal) <- 1

    C_epistatic <- matrix(epistatic_correlation, ncol = d, nrow = d)
    diag(C_epistatic) <- 1

    C_error <- matrix(0, ncol = d, nrow = d)
    diag(C_error) <- 1

    log$debug("Desired marginal correlation: %f", marginal_correlation)
    beta <- mvtnorm::rmvnorm(n_causal, sigma = C_marginal)
    log$debug("Correlation of simulated marginal effects: %s", cor(beta))

    n_epistatic_effects <- n_group1_trait * (n_trait_specific - n_group1_trait) +
        ncol(X_epi_pleio)
    log$debug("Number of epistatic effects: %s", n_epistatic_effects)
    if (n_epistatic_effects > 0) {
        log$debug("Desired epistatic correlation: %f", epistatic_correlation)
        alpha <- mvtnorm::rmvnorm(n_epistatic_effects, sigma = C_epistatic)
        log$debug("Correlation of simulated epistatic effects: %s", cor(alpha))
    } else {
        log$debug("Return empty effect matrix.")
        alpha <- matrix(0, ncol = d, nrow = 0)
    }
    log$debug("Desired error correlation: %f", 0)
    error <- mvtnorm::rmvnorm(n_samples, sigma = C_error)
    log$debug("Correlation of simulated error: %s", cor(error))
    snp.ids.trait <- setdiff(snp.ids.filtered, pleiotropic_set)

    pleiotropic_interactions <- tibble(
        group1 = rep(pleio_split$group1, each = length(pleio_split$group2)),
        group2 = rep(pleio_split$group2, times = length(pleio_split$group1))
    )

    interactions <- tibble()
    additive <- tibble()

    for (j in 1:d) {
        ## select causal SNPs
        log$debug("Simulating phenotype %d", j)
        causal_snps_j <- sample(snp.ids.trait, n_causal - n_pleiotropic, replace = F)
        trait_specific_additive <- c(causal_snps_j, pleiotropic_set)
        trait_specific_snps <- sample(causal_snps_j, n_trait_specific, replace = F)
        trait_grouped <- split(trait_specific_snps, f_trait)
        trait_specific_j_1 <- trait_grouped$group1
        trait_specific_j_2 <- trait_grouped$group2
        log$debug("Length causal set: %d", length(causal_snps_j))
        log$debug("Length pleiotropic set: %d", length(pleiotropic_set))
        log$debug("Length trait specific set 1: %d", length(trait_specific_j_1))
        log$debug("Length trait specific set 2: %d", length(trait_specific_j_2))

        log$debug("Head of causal SNPs: %s", head(trait_specific_additive))
        log$debug("Head of trait specific SNPs group 1: %s", head(trait_specific_j_1))
        log$debug("Head of trait specific SNPs group 2: %s", head(trait_specific_j_2))

        # create trait_specific interaction matrix
        X_causal_j <- X[, c(causal_snps_j, pleiotropic_set)]  # all SNPs have additive effects
        X_trait_specific_j_1 <- X[, trait_specific_j_1]
        X_trait_specific_j_2 <- X[, trait_specific_j_2]

        trait_specific_interactions <- tibble(
            group1 = rep(trait_specific_j_1, each = length(trait_specific_j_2)),
            group2 = rep(trait_specific_j_2, times = length(trait_specific_j_1))
        )

        start_interactions <- proc.time()
        log$debug("Computing interactions. This may take a while.")
        X_epi <- foreach(i = seq_len(length(trait_specific_j_1)), .combine = cbind) %do% {
                X_trait_specific_j_1[, i] * X_trait_specific_j_2
        }
        X_epi <- cbind(X_epi_pleio, X_epi)

        time_interactions <- proc.time() - start_interactions
        log$debug("Interactions X_epi computed in %f", time_interactions[3])
        log$debug(
            "Dimension of interaction matrix X_epi: %d x %d", nrow(X_epi),
            ncol(X_epi)
        )

        # marginal effects
        X_marginal <- X_causal_j
        beta_j <- beta[, j]
        y_marginal <- X_marginal %*% beta_j
        y_marginal <- y_marginal * sqrt(heritability[j] * rho/c(var(y_marginal)))
        log$debug("Variance scaled y_marginal: %f", var(y_marginal))

        # pairwise epistatic effects
        alpha_j <- alpha[, j]
        if (n_epistatic_effects > 0) {
            y_epi <- X_epi %*% alpha_j
            y_epi <- y_epi * sqrt(heritability[j] * (1 - rho)/c(var(y_epi)))
            log$debug("Variance scaled y_epi: %f", var(y_epi))
            trait_interactions <- bind_rows(pleiotropic_interactions,
                                            trait_specific_interactions) %>%
                mutate(effect_size = alpha_j) %>%
                mutate(trait = j)
            interactions <- bind_rows(interactions, trait_interactions)
        } else {
            y_epi <- 0 * y_marginal
            log$debug("y_epi vector of zeros size y_marginal: %s", y_epi)
        }
        trait_additive <- tibble(
            id = trait_specific_additive,
            effect_size = beta_j,
            trait = j
        )
        additive <- bind_rows(additive, trait_additive)

        # unexplained phenotypic variation
        y_err <- error[, j]
        y_err <- y_err * sqrt((1 - heritability[j])/c(var(y_err)))

        y <- y_marginal + y_epi + y_err
        Y <- cbind(Y, y)
        causal_snps[[paste0("phenotype_", j)]] <- list(
            causal_snps = c(causal_snps_j, pleiotropic_set),
            pleiotropic_groups = pleio_split, trait_specific_groups = trait_grouped,
            alpha = alpha_j, beta = beta_j
        )
    }

    if (n_epistatic_effects > 0) {
        epistatic <- tidyr::pivot_longer(
                        interactions,
                        cols = c("group1", "group2"),
                        names_to = "group",
                        names_prefix = "group",
                        values_to = "id"
                        ) %>%
                     mutate(name = sprintf("snp_%05d", id)) %>%
                     dplyr::group_by(trait, id, name) %>%
                     summarise(total_effect = sum(abs(effect_size))) %>%
                     mutate(pleiotropic = (id  %in% pleiotropic_set))
    } else {
        epistatic <- NULL
        interactions <- NULL
    }

    colnames(genotype_matrix) <- seq_len(ncol(genotype_matrix)) %>%
        sprintf(fmt = "snp_%05d")  # column names names for SNPs
    colnames(Y) <- seq_len(ncol(Y)) %>%
        sprintf(fmt = "p_%02d")  # column names names for phenotypes

    log$debug("Phenotype data: %s", head(Y))
    log$debug("Phenotype correlation: %s", cor(Y))
    log$debug("Correlation of simulated epistatic effects: %s", cor(alpha))

    # return data
    parameter_names <- c("number_samples",
                         "number_snps",
                         "number_phenotypes",
                         "number_causal_snps",
                         "number_epistatic_effects",
                         "number_pleiotropic_snps",
                         "number_trait_specific_snps",
                         "heritability",
                         "rho",
                         "epistatic_correlation",
                         "marginal_correlation",
                         "group_ratio_pleiotropic",
                         "group_ratio_trait",
                         "seed"
    )
    parameters <- tibble()
    for (j in seq_len(d)) {
        parameter_values <- c(n_samples,
                              n_snp,
                              d,
                              n_causal,
                              n_epistatic_effects,
                              n_pleiotropic,
                              n_trait_specific,
                              heritability[j],
                              rho,
                              epistatic_correlation,
                              marginal_correlation,
                              group_ratio_pleiotropic,
                              group_ratio_trait,
                              seed)
        parameters <- bind_rows(parameters,
                                tidyr::tibble(name = parameter_names,
                                    value = parameter_values)
                                )
    }
    parameters <- parameters %>%
        mutate(trait = rep(colnames(Y), each = length(parameter_names)))
    simulated_pleiotropic_epistasis_data <- list(
        parameters = parameters, phenotype = Y, genotype = genotype_matrix,
        additive = additive, epistatic = epistatic, interactions = interactions,
        snps.filtered = snp.ids.filtered
    )

    return(simulated_pleiotropic_epistasis_data)
}

get_factors <- function(n1, n) {
    return(
        c(
            rep("group1", n1),
            rep("group2", n - n1)
        )
    )
}

