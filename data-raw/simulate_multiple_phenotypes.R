set.seed(982348)
library('MASS')

file_name <- 'data/random_genotype_matrix.rda'
load(file_name)

X <- scale(random_genotype_matrix)

n_samples <- nrow(X) # number of genotype samples
n_snp <- ncol(X) # number of SNPs
H2 <- 0.6 # pve: phenotypic variance explained/broad-sense heritability (H^2)
rho <- 0.5 # rho: measures the portion of H^2 that is contributed by the marginal (additive) effects
d <- 5 # number of phenotypes
marginal_correlation <- 0.3
epistatic_correlation <- 0.3

n_causal <- 20 # number of SNPs to be causal in every phenotype
n_causal_epi <- 3 # number of epistatic causal SNPs slected per interaction group and phenotype
n_causal_pleio <- 5 # number of SNPs to be involved in pleiotropic effects in every phenotype

snp.ids <- 1:n_snp

Y_marginal <- c()
Y_epistatic <- c()
Y_error <- c()
causal_snps <- list()
pleiotropic_set <- sample(snp.ids, n_causal_pleio, replace = F) # declare peleiotropic SNPs before since they have to be present in every phenotype

for (j in 1:d) {
  ## select causal SNPs
  causal_snps_j <- sample(snp.ids[-pleiotropic_set], n_causal, replace = F)
  epistatic_set_j_1 <- sample(causal_snps_j, n_causal_epi, replace = F)
  epistatic_set_j_2 <- sample(causal_snps_j[! causal_snps_j %in% epistatic_set_j_1], n_causal_epi, replace = F)

  epistatic_set_j_1 <- c(epistatic_set_j_1, pleiotropic_set) # the epistatic pleiotropic effects are included in epistatic interaction group 1

  # create epistatic interaction matrix
  X_causal_j <- X[, c(causal_snps_j, pleiotropic_set)] # all SNPs have additive effects
  X_epistatic_j_1 <- X[, epistatic_set_j_1]
  X_epistatic_j_2 <- X[, epistatic_set_j_2]
  X_epi <- c()
  for (i in 1:length(epistatic_set_j_1)) {
    X_epi <- cbind(X_epi, X_epistatic_j_1[, i] * X_epistatic_j_2)
  }
  dim(X_epi)

  # marginal effects
  X_marginal <- X_causal_j
  beta <- rnorm(dim(X_marginal)[2])
  y_marginal <- X_marginal %*% beta
  beta = beta * sqrt(H2 * rho / c(var(y_marginal)))
  y_marginal=X_marginal %*% beta

  # pairwise epistatic effects
  alpha <- rnorm(dim(X_epi)[2])
  y_epi <- X_epi %*% alpha
  alpha = alpha * sqrt(H2 * (1 - rho) / c(var(y_epi)))
  y_epi = X_epi %*% alpha

  # unexplained phenotypic variation
  y_err <- rnorm(n_samples)
  y_err <- y_err * sqrt((1 - H2) / c(var(y_err)))

  Y_marginal <- cbind(Y_marginal, y_marginal)
  Y_epistatic <- cbind(Y_epistatic, y_epi)
  Y_error <- cbind(Y_error, y_err)
  causal_snps[[paste0('phenotype_', j)]] <- list(
    'causal_snps' = causal_snps_j,
    'epistatic_1' = epistatic_set_j_1,
    'epistatic_2' = epistatic_set_j_2,
    'alpha' = alpha,
    'beta' = beta
  )
}

# scale marginal data variance and correlation
C <- matrix(marginal_correlation, ncol = d, nrow = d)
diag(C) <- 1
Y_marginal <- Y_marginal %*% chol(C)

# scale epistatic data variance and correlation
C <- matrix(epistatic_correlation, ncol = d, nrow = d)
diag(C) <- 1
Y_epistatic <- Y_epistatic %*% chol(C)

# scale error variance
C <- matrix(0, ncol = d, nrow = d)
diag(C) <- 1
Y_error <- Y_error %*% chol(C)

Y <- Y_marginal + Y_epistatic + Y_error

colnames(X) <- seq_len(ncol(X)) %>% sprintf(fmt = "snp_%05d") # column names names for SNPs
colnames(Y) <- seq_len(ncol(Y)) %>% sprintf(fmt = "p_%02d") # column names names for phenotypes

# export data
simulated_pleiotropic_epistasis_data <- list(
  number_snp = n_snp,
  number_samples = n_samples,
  pve = H2,
  rho = rho,
  phenotype = Y,
  genotype = file_name,
  causal_snps = causal_snps
)

usethis::use_data(simulated_pleiotropic_epistasis_data, overwrite = TRUE)
