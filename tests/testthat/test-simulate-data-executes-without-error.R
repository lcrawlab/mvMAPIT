test_that("Simulate multiple phenotypes returns apropriate phenotype object", {
  # given
  p <- 10
  n <- 5
  d <- 3
  X <- matrix(runif(p * n), ncol = p)
  
  # when
  data <- simulate_phenotypes(X, d = d, logLevel = 'DEBUG')
  
  # then
  expect_equal(nrow(data$phenotype), n)
  expect_equal(ncol(data$phenotype), d)
  expect_equal(sum(is.na(data$phenotype)), 0) # no NA values in phenotype
  expect_equal(length(data), 7)
})

test_that("Simulate multiple phenotypes returns causal SNPs", {
  # given
  p <- 10
  f <- 0.4
  g <- 0.1
  n <- 5
  d <- 3
  X <- matrix(runif(p * n), ncol = p)
  total_causal <- (f * p) + ceiling((f * p) * g) # single trait SNPs plus pleiotropic SNPs
  
  # when
  data <- simulate_phenotypes(X, causal_fraction = f, pleiotropic_fraction = g, d = d, logLevel = 'DEBUG')
  
  # then
  expect_equal(length(data$causal_snps$phenotype_1), total_causal)
})
