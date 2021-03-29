library(dplyr)

# As in the paper, we randomly choose causal variants that classify into three groups:
#
#   (i) a small set of interaction SNPs,
#   (ii) a larger set of interaction SNPs, and
#   (iii) a large set of additive SNPs.
#
# In the simulations carried out in this study, SNPs interact between sets,
# so that SNPs in the first group interact with SNPs in the second group,
# but do not interact with variants in their own group (the same applies to the
# second group). One may view the SNPs in the first set as the “hubs” in an
# interaction map. We are reminded that interaction (epistatic) effects are
# different from additive effects. All causal SNPs in both the first and second
# groups have additive effects and are involved in pairwise interactions, while 
# causal SNPs in the third set only have additive effects.

set.seed(11151990)

n_samples <- 3e2
n_snp <- 1e3
n_causal_snp <- 1e2 # has to be less than n_snp

n_causal_1 <- 10 # Set 1 of causal SNPs
n_causal_2 <- 10 # Set 2 of Causal SNPs
n_causal_3 <- n_causal_snp - n_causal_1 - n_causal_2 # Set 3 of Causal SNPs with only marginal effects

pve <- 0.6 # pve: phenotypic variance explained/broad-sense heritability (H^2)
rho <- 0.5 # rho: measures the portion of H^2 that is contributed by the marginal (additive) effects

# simulate data with minor allele frequency > 0.05
maf <- 0.05 + 0.45 * runif(n_samples * n_snp)
genotypes <- (runif(n_samples * n_snp) < maf) + (runif(n_samples * n_snp) < maf) # simulate two independent haplotype vectors and sum them for genotype
genotype_matrix <- matrix(
  as.double(genotypes),
  nrow = n_samples,
  ncol = n_snp,
  byrow = TRUE
)
X_mean <- apply(genotype_matrix, 2, mean) # mean over SNPs
X_sd <- apply(genotype_matrix, 2, sd)
X <- t((t(genotype_matrix) - X_mean) / X_sd) # X: centered and scaled genotype matrix

n_snp
snp.ids <- 1:n_snp
# select causal SNPs
id_causal_1 <- sample(snp.ids, n_causal_1, replace = F)
id_causal_2 <- sample(snp.ids[-id_causal_1], n_causal_2, replace = F)
id_causal_3 <- sample(snp.ids[-c(id_causal_1, id_causal_2)], n_causal_3, replace = F)

X_causal1 <- X[, id_causal_1]
X_causal2 <- X[, id_causal_2]
X_causal3 <- X[, id_causal_3]

# pairwise interactions of causal SNPs from group 1 with causaul SNPs from group 2
W <- c()
for (i in 1:n_causal_1) {
  # element wise multiplication of all sample genotypes at SNP i with all sample genotypes of SNPs in group 2
  W <- cbind(W, X_causal1[, i] * X_causal2)
}

X_marginal <- cbind(X_causal1, X_causal2, X_causal3)
beta <- rnorm(dim(X_marginal)[2]) # random effect size for causal SNPs additive effect
y_marginal <- X_marginal %*% beta # phenotype due to linear effects of causal SNPs unscaled
beta <- beta * sqrt(pve * rho / c(var(y_marginal))) # rescale effect to the desired proportion of variance explained due to linear effects
y_marginal <- X_marginal %*% beta # phenotype due to linear effects of causal SNPs

# repeat for epistatic effect sizes
alpha <- rnorm(dim(W)[2])
y_epi <- W %*% alpha
alpha <- alpha * sqrt(pve * (1 - rho) / c(var(y_epi)))
y_epi <- W %*% alpha

# remaining variance dure to random error
y_err <- rnorm(n_samples)
y_err <- y_err * sqrt((1 - pve) / c(var(y_err)))

y <- y_marginal + y_epi + y_err # phenotype with desired proportions of variance


colnames(X) <- seq_len(ncol(X)) %>% sprintf(fmt = "snp_%05d") # column names names for SNPs

epistatic_snps <- colnames(X)[c(id_causal_1, id_causal_2)]
additive_snps <- colnames(X)[id_causal_3]

# export data
simulated_epistasis_data <- list(
  number_snp = n_snp,
  number_samples = n_samples,
  number_causal_snp = n_causal_snp,
  pve = pve,
  rho = rho,
  phenotype = y,
  genotype = X,
  linear_effects = beta,
  interaction_effects = alpha,
  epistatic_snps = epistatic_snps,
  additive_snps = additive_snps
)

usethis::use_data(simulated_epistasis_data, overwrite = TRUE)
