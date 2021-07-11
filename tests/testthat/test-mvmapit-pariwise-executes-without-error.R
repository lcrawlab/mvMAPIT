test_that("hybrid = FALSE. phenotypeCovariance = ''", {
  # given
  p <- 2
  n <- 10
  d <- 3
  pvalues <- matrix(c(0.4705, 0.5025, 0.4007, 0.3927, 0.6591, 0.6629, 0.4171, 0.0285, 0.3101, 0.1738, 0.0308, 0.0503),
                    nrow = p, ncol = 6)
  set.seed(853)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   hybrid = FALSE,
                   cores = 1,
                   phenotypeCovariance = '',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-4)
})

test_that("hybrid = FALSE, test = davies. phenotypeCovariance = ''", {
  # given
  p <- 2
  n <- 10
  d <- 3
  pvalues <- matrix(
                    c(0.01559003,
                      0.01606844,
                      0.0,
                      0.0,
                      0.3795864,
                      0.3622176,
                      0.0002069482,
                      0.0962537846,
                      4.045376e-07,
                      1.987534e-05,
                      0.3353284,
                      0.2971684),
                    nrow = p, ncol = 6)
  set.seed(853)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   hybrid = FALSE,
                   test = 'davies',
                   cores = 1,
                   phenotypeCovariance = '',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-4)
})


test_that("hybrid = TRUE. phenotypeCovariance = ''", {
  # given
  p <- 2
  n <- 10
  d <- 3
  pvalues <- matrix(
                    c(0.4705429,
                      0.5024536,
                      0.4007342,
                      0.3926739,
                      0.6591476,
                      0.6629176,
                      0.417081,
                      0.0963,
                      0.3101424,
                      0.1737842,
                      0.3353,
                      0.05034315),
                    nrow = p, ncol = 6)
  set.seed(853)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   hybrid = TRUE,
                   cores = 1,
                   phenotypeCovariance = '',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-4)
})