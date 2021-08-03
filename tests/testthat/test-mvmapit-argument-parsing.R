test_that("MvMapit can take a vector as phenotype input. hybrid = FALSE, test = normal", {
  # given
  p <- 10
  n <- 4
  pvalues <- matrix(rep(0.48001, 10), ncol = 1)
  set.seed(853)
  X <- matrix(runif(p * n), ncol = p)
  Y <- c(runif(n))
  # when
  mapit <- MvMAPIT(t(X),
                   Y,
                   accuracy = 1e-5,
                   cores = 1,
                   phenotypeCovariance = 'identity',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-3)
})

test_that("MvMapit can take a vector as phenotype input. hybrid = FALSE, test = davies", {
  # given
  p <- 10
  n <- 4
  pvalues <- matrix(rep(0.209, 10), ncol = 1)
  set.seed(853)
  X <- matrix(runif(p * n), ncol = p)
  Y <- c(runif(n))
  # when
  mapit <- MvMAPIT(t(X),
                   Y,
                   test = 'davies',
                   accuracy = 1e-5,
                   cores = 1,
                   phenotypeCovariance = 'identity',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-3)
})
