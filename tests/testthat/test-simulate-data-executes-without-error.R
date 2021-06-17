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
  expect_equal(length(data), 8)
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
