test_that("Simulate multiple phenotypes returns apropriate phenotype object", {
  # given
  p <- 20
  f <- 10
  g <- 4
  n <- 5
  d <- 3
  X <- matrix(runif(p * n), ncol = p)

  # when
  data <- simulate_phenotypes(X,
                              n_causal = f,
                              n_trait_specific = g,
                              n_pleiotropic = g,
                              d = d,
                              maf_threshold = 0.0,
                              logLevel = 'ERROR')

  # then
  expect_equal(nrow(data$phenotype), n)
  expect_equal(ncol(data$phenotype), d)
  expect_equal(sum(is.na(data$phenotype)), 0) # no NA values in phenotype
  expect_equal(length(data), 8)
})

test_that("Simulate multiple phenotypes returns causal SNPs", {
  # given
  p <- 20
  f <- 10
  g <- 4
  n <- 5
  d <- 3
  X <- matrix(runif(p * n), ncol = p)
  total_causal <- (f) # single trait SNPs plus pleiotropic SNPs

  # when
  data <- simulate_phenotypes(X,
                              n_causal = f,
                              n_trait_specific = g,
                              n_pleiotropic = g,
                              d = d,
                              maf_threshold = 0.0,
                              logLevel = 'ERROR')

  # then
  expect_equal(length(data$snps$phenotype_1$causal_snps), total_causal)
})

test_that("Simulate multiple phenotypes remove SNPs with low maf", {
  # given
  p <- 20
  f <- 10
  g <- 4
  n <- 5
  d <- 3
  maf <- 0.05 + 0.45 * runif(p)
  set.seed(12345)
  X <- matrix((runif(p * n) >= maf) + (runif(p * n) >= maf), ncol = p)
  X[, 1] <- 0
  X[, 13] <- 0
  total_causal <- (f * p) # single trait SNPs plus pleiotropic SNPs

  # when
  data <- simulate_phenotypes(X,
                              n_causal = f,
                              n_trait_specific = g,
                              n_pleiotropic = g,
                              d = d,
                              maf_threshold = 0.05,
                              logLevel = 'ERROR')
  data$snps.filtered

  # then
  expect_true(!(1 %in% data$snps.filtered) & !(13 %in% data$snps.filtered))
})

test_that("Simulation of groups with given ratio works as expected.", {
  # given
  p <- 50
  f <- 30
  g <- 8
  h <- 10
  n <- 5
  d <- 3
  X <- matrix(runif(p * n), ncol = p)
  total_causal <- (f) # single trait SNPs plus pleiotropic SNPs
  correct_1 <- 2
  correct_2 <- 6
  correct_3 <- 2
  correct_4 <- 8

  # when
  data <- simulate_phenotypes(X,
                              n_causal = f,
                              n_trait_specific = g,
                              n_pleiotropic = h,
                              d = d,
                              group_ratio_trait = 3,
                              group_ratio_pleiotropic = 4,
                              maf_threshold = 0.0,
                              logLevel = 'ERROR')

  # then
  expect_equal(length(data$snps$phenotype_1$trait_specific_groups$group1), correct_1)
  expect_equal(length(data$snps$phenotype_1$trait_specific_groups$group2), correct_2)
  expect_equal(length(data$snps$phenotype_1$pleiotropic_groups$group1), correct_3)
  expect_equal(length(data$snps$phenotype_1$pleiotropic_groups$group2), correct_4)
})

test_that("simulate_phenotypes can handle zero size bn_trait_specific groups.", {
  # given
  p <- 50
  f <- 30
  g <- 0
  h <- 10
  n <- 5
  d <- 3
  X <- matrix(runif(p * n), ncol = p)
  total_causal <- (f) # single trait SNPs plus pleiotropic SNPs
  correct_1 <- 0
  correct_2 <- 0

  # when
  data <- simulate_phenotypes(X,
                              n_causal = f,
                              n_trait_specific = g,
                              n_pleiotropic = h,
                              d = d,
                              group_ratio_trait = 3,
                              group_ratio_pleiotropic = 4,
                              maf_threshold = 0.0,
                              logLevel = 'ERROR')

  # then
  expect_equal(length(data$snps$phenotype_1$trait_specific_groups$group1), correct_1)
  expect_equal(length(data$snps$phenotype_1$trait_specific_groups$group2), correct_1)
})

test_that("simulate_phenotypes can handle zero size bn_pleiotropic groups.", {
  # given
  p <- 50
  f <- 30
  g <- 10
  h <- 0
  n <- 5
  d <- 3
  X <- matrix(runif(p * n), ncol = p)
  total_causal <- (f) # single trait SNPs plus pleiotropic SNPs
  correct_1 <- 0
  correct_2 <- 0

  # when
  data <- simulate_phenotypes(X,
                              n_causal = f,
                              n_trait_specific = g,
                              n_pleiotropic = h,
                              d = d,
                              group_ratio_trait = 3,
                              group_ratio_pleiotropic = 4,
                              maf_threshold = 0.0,
                              logLevel = 'DEBUG')

  # then
  expect_equal(length(data$snps$phenotype_1$pleiotropic_groups$group1), correct_1)
  expect_equal(length(data$snps$phenotype_1$pleiotropic_groups$group2), correct_2)
})

test_that("test run", {
ind <- 1e2
nsnp <- 100
H2 <- 0.6
rho <- 0.5
maf <- 0.05 + 0.45*runif(nsnp)
X   <- (runif(ind*nsnp) < maf) + (runif(ind*nsnp) < maf)
X   <- matrix(as.double(X),ind,nsnp,byrow = TRUE)
s <- 95345 # sample.int(10000, 1)
sim <- simulate_phenotypes(X,
                           n_causal = 30,
                           n_pleiotropic = 6,
                           n_trait_specific = 4,
                           epistatic_correlation = 0.8,
                           H2 = H2, rho = rho, logLevel = 'ERROR', seed = s)
})
