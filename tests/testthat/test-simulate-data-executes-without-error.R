test_that("Simulate multiple phenotypes returns apropriate phenotype object", {
  # given
  p <- 20
  f <- 0.4
  g <- 0.5
  n <- 5
  d <- 3
  X <- matrix(runif(p * n), ncol = p)

  # when
  data <- simulate_phenotypes(X, causal_fraction = f, pleiotropic_fraction = g, d = d, maf_threshold = 0.0, logLevel = 'ERROR')

  # then
  expect_equal(nrow(data$phenotype), n)
  expect_equal(ncol(data$phenotype), d)
  expect_equal(sum(is.na(data$phenotype)), 0) # no NA values in phenotype
  expect_equal(length(data), 9)
})

test_that("Simulate multiple phenotypes returns causal SNPs", {
  # given
  p <- 20
  f <- 0.4
  g <- 0.5
  n <- 5
  d <- 3
  X <- matrix(runif(p * n), ncol = p)
  total_causal <- (f * p) # single trait SNPs plus pleiotropic SNPs

  # when
  data <- simulate_phenotypes(X, causal_fraction = f, pleiotropic_fraction = g, d = d, maf_threshold = 0.0, logLevel = 'ERROR')

  # then
  expect_equal(length(data$snps$phenotype_1$causal_snps), total_causal)
})

test_that("Simulate multiple phenotypes remove SNPs with low maf", {
  # given
  p <- 20
  f <- 0.4
  g <- 0.5
  n <- 5
  d <- 3
  maf <- 0.05 + 0.45 * runif(p)
  set.seed(12345)
  X <- matrix((runif(p * n) >= maf) + (runif(p * n) >= maf), ncol = p)
  X[, 1] <- 0
  X[, 13] <- 0
  total_causal <- (f * p) # single trait SNPs plus pleiotropic SNPs

  # when
  data <- simulate_phenotypes(X, causal_fraction = f, pleiotropic_fraction = g, d = d, maf_threshold = 0.05, logLevel = 'ERROR')
  data$snps.filtered

  # then
  expect_true(!(1 %in% data$snps.filtered) & !(13 %in% data$snps.filtered))
})
