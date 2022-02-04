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
#' @param causal_fraction Fraction of the SNPs that are causal.
#' @param epistatic_fraction Fraction of the causal SNPs with single trait epistatic effects.
#' @param pleiotropic_fraction Fraction of SNPs with pleiotropic effects.
#' @param H2 Broad-sense heritability.
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
simulate_phenotypes <- function(genotype_matrix,
                                causal_fraction = 0.2,
                                epistatic_fraction = 0.1,
                                pleiotropic_fraction = 0.1,
                                H2 = 0.6,
                                d = 2,
                                rho = 0.8,
                                marginal_correlation = 0.3,
                                epistatic_correlation = 0.3,
                                maf_threshold = 0.01,
                                seed = 67132,
                                logLevel = 'INFO',
                                logFile = NULL) {
  set.seed(seed)
  coll <- makeAssertCollection()
  assertDouble(causal_fraction, lower = 0, upper = 1, add = coll)
  assertDouble(epistatic_fraction, lower = 0, upper = 0.3, add = coll)
  assertDouble(pleiotropic_fraction, lower = 0, upper = 0.5, add = coll)
  assertDouble(H2, lower = 0, upper = 1, add = coll)
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
  log <- logging::getLogger('simulate_phenotypes')
  if(!is.null(logFile)) {
    filePath <- file.path(getwd(),logFile)
    log$debug('Logging to file: %s', filePath)
    log$addHandler(logging::writeToFile, file=filePath)
  }

  snp.ids <- 1:ncol(genotype_matrix)
  maf <- colMeans(genotype_matrix) / 2
  X <- genotype_matrix
  maf_compliant <- (maf > maf_threshold) & (maf < 1 - maf_threshold)
  # scale produces NaN when the columns have zero variance
  snp.ids.filtered <- snp.ids[complete.cases(t(scale(X))) & maf_compliant]

  n_samples <- nrow(X) # number of genotype samples
  n_snp <- ncol(X) # number of SNPs passing quality control
  log$debug('Scaled genotype matrix: %d x %d', n_samples, n_snp)
  log$debug('Disregard %d variants due to zero variance or small minor allele frequency.', ncol(genotype_matrix) - length(snp.ids.filtered))
  log$debug('Minor allele frequency threshold %f.', maf_threshold)

  n_causal <- ceiling(n_snp * causal_fraction) # number of SNPs to be causal in every phenotype
  n_causal_epi <- ceiling(n_causal * epistatic_fraction) # number of epistatic causal SNPs slected per interaction group and phenotype
  n_causal_pleio <- ceiling(n_causal * pleiotropic_fraction) # number of SNPs to be involved in pleiotropic effects in every phenotype
  log$debug('Number of causal SNPs: %d', n_causal)
  log$debug('Number of epistatic SNPs: %d', n_causal_epi)
  log$debug('Number of pleiotropic SNPs: %d', n_causal_pleio)
  log$debug('NA in raw genotype matrix: %d', sum(is.na(genotype_matrix)))
  log$debug('NA in scaled genotype matrix: %d', sum(is.na(X)))

  Y <- c()
  causal_snps <- list()

  pleiotropic_set <- sample(snp.ids.filtered, n_causal_pleio, replace = F) # declare peleiotropic SNPs before since they have to be present in every phenotype
  pleio_split <- split(pleiotropic_set, f = c('a', 'b'))
  X_pleio_a <- X[, pleio_split$a]
  X_pleio_b <- X[, pleio_split$b]
  X_epi_pleio <- foreach(i=seq_len(length(pleio_split$a)), .combine=cbind) %do% {
    X_pleio_a[, i] * X_pleio_b
  }
  log$debug('Dimensions of pleiotropic interaction matrix: %d x %d', nrow(X_epi_pleio), ncol(X_epi_pleio))

  log$debug('Draw effects from multivariate normal with desired correlation.')

  C_marginal <- matrix(marginal_correlation, ncol = d, nrow = d)
  diag(C_marginal) <- 1

  C_epistatic <- matrix(epistatic_correlation, ncol = d, nrow = d)
  diag(C_epistatic) <- 1

  C_error <- matrix(0, ncol = d, nrow = d)
  diag(C_error) <- 1

  log$debug('Desired marginal correlation: %f', marginal_correlation)
  beta <- mvtnorm::rmvnorm(n_causal, sigma = C_marginal)
  log$debug('Correlation of simulated marginal effects: %s', cor(beta))

  log$debug('Desired epistatic correlation: %f', epistatic_correlation)
  alpha <- mvtnorm::rmvnorm(n_causal_epi * n_causal_pleio + ncol(X_epi_pleio),
                            sigma = C_epistatic)
  log$debug('Correlation of simulated epistatic effects: %s', cor(alpha))

  log$debug('Desired error correlation: %f', 0)
  error <- mvtnorm::rmvnorm(n_samples, sigma = C_error)
  log$debug('Correlation of simulated error: %s', cor(error))

  for (j in 1:d) {
    ## select causal SNPs
    log$debug('Simulating phenotype %d', j)
    causal_snps_j <- sample(snp.ids.filtered[-pleiotropic_set], n_causal - n_causal_pleio, replace = F)
    epistatic_set_j_1 <- pleiotropic_set # the epistatic pleiotropic effects are included in epistatic interaction group 1
    epistatic_set_j_2 <- sample(causal_snps_j, n_causal_epi, replace = F)

    log$debug('Head of causal SNPs: %s', head(c(causal_snps_j, pleiotropic_set)))
    log$debug('Head of epistatic SNPs group 1: %s', head(epistatic_set_j_1))
    log$debug('Head of epistatic SNPs group 2: %s', head(epistatic_set_j_2))

    # create epistatic interaction matrix
    X_causal_j <- X[, c(causal_snps_j, pleiotropic_set)] # all SNPs have additive effects
    X_epistatic_j_1 <- as.matrix(X[, epistatic_set_j_1])
    X_epistatic_j_2 <- X[, epistatic_set_j_2]

    start_interactions <- proc.time()
    log$debug('Computing interactions. This may take a while.')
    X_epi <- foreach(i=seq_len(length(epistatic_set_j_1)), .combine=cbind) %do% {
      X_epistatic_j_1[, i] * X_epistatic_j_2
    }
    X_epi <- cbind(X_epi_pleio, X_epi)

    time_interactions <- proc.time() - start_interactions
    log$debug('Interactions X_epi computed in %f', time_interactions[3])
    log$debug('Dimension of interaction matrix X_epi: %d x %d', nrow(X_epi), ncol(X_epi))

    # marginal effects
    X_marginal <- X_causal_j
    beta_j <- beta[, j]
    y_marginal <- X_marginal %*% beta_j
    y_marginal <- y_marginal * sqrt(H2 * rho / c(var(y_marginal)))
    log$debug('Variance scaled y: %f', var(y_marginal))

    # pairwise epistatic effects
    alpha_j <- alpha[, j]
    y_epi <- X_epi %*% alpha_j
    y_epi <- y_epi * sqrt(H2 * (1 - rho) / c(var(y_epi)))
    log$debug('Variance scaled y_epi: %f', var(y_epi))

    # unexplained phenotypic variation
    y_err <- error[, j]
    y_err <- y_err * sqrt((1 - H2) / c(var(y_err)))

    y <- y_marginal + y_epi + y_err
    Y <- cbind(Y, y)
    causal_snps[[paste0('phenotype_', j)]] <- list(
      'causal_snps' = c(causal_snps_j, pleiotropic_set),
      'epistatic_1' = epistatic_set_j_1,
      'epistatic_2' = epistatic_set_j_2,
      'pleiotropic' = pleiotropic_set,
      'alpha' = alpha_j,
      'beta' = beta_j
    )
  }


  colnames(genotype_matrix) <- seq_len(ncol(genotype_matrix)) %>% sprintf(fmt = "snp_%05d") # column names names for SNPs
  colnames(Y) <- seq_len(ncol(Y)) %>% sprintf(fmt = "p_%02d") # column names names for phenotypes

  log$debug('Phenotype data: %s', head(Y))
  log$debug('Phenotype correlation: %s', cor(Y))
  log$debug('Correlation of simulated epistatic effects: %s', cor(alpha))

  # return data
  simulated_pleiotropic_epistasis_data <- list(
    number_samples = n_samples,
    pve = H2,
    rho = rho,
    phenotype = Y,
    genotype = genotype_matrix,
    snps = causal_snps,
    snps.filtered = snp.ids.filtered,
    seed = seed
  )

  return(simulated_pleiotropic_epistasis_data)
}
