test_that("MvMapit executes without error when hybrid = TRUE.", {
  # given
  p <- 10
  n <- 5
  d <- 1
  pvalues <- c(0.6772546,
               0.5690567,
               0.3366181,
               0.5523170,
               0.6307954,
               0.4436580,
               0.7090439,
               0.6147064,
               0.7128274,
               0.6088642)
  set.seed(5)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   cores = 1,
                   phenotypeCovariance = 'covariance',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-7)
})

test_that("MvMapit executes without error when hybrid = FALSE.", {
  # given
  p <- 10
  n <- 5
  d <- 3
  pvalues <- c(0.3409648,
               0.5879707,
               0.5973206,
               0.2912437,
               0.9626492,
               0.9496355,
               0.5612714,
               0.9947883,
               0.6359662,
               0.8031088)
  set.seed(6)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X), 
                   t(Y), 
                   hybrid = FALSE, 
                   cores = 1, 
                   phenotypeCovariance = 'covariance',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-7)
})

test_that("MvMapit executes without error when hybrid = FALSE and test = davies.", {
  # given
  p <- 10
  n <- 5
  d <- 3
  pvalues <- c(0.4855836,
               0.5042975,
               0.6235587,
               0.2560563,
               0.9150879,
               0.8862566,
               0.6406923,
               0.9294568,
               0.5974718,
               0.9704179)
  set.seed(6)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X), 
                   t(Y), 
                   hybrid = FALSE, 
                   test = 'davies', 
                   cores = 1, 
                   variantIndex = c(1:p),
                   phenotypeCovariance = 'covariance',
                   logLevel = "ERROR")
  # then
  expect_equal(length(mapit$pvalues), p)
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-7)
})

test_that("MvMapit times computations when hybrid = TRUE.", {
  # given
  p <- 10
  n <- 5
  d <- 1
  set.seed(5)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   cores = 1,
                   hybrid = TRUE,
                   phenotypeCovariance = 'covariance',
                   logLevel = "ERROR")
  # then
  expect_equal(length(mapit$timings), 6)
})

test_that("MvMapit executes without error when C is not NULL, hybrid = TRUE.", {
  # given
  p <- 10
  n <- 5
  d <- 1
  pvalues <- c(0.6828469,
               0.6014190,
               0.7815753,
               0.8516215,
               0.5315767,
               0.9810549,
               0.4784025,
               0.3926182,
               0.7893141,
               0.9439612)
  set.seed(6)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  C <- matrix(runif(n * n), ncol = n)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   C = C,
                   cores = 1,
                   hybrid = TRUE,
                   phenotypeCovariance = 'covariance',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-7)
})


test_that("hybrid = FALSE. phenotypeCovariance = identity", {
  # given
  p <- 10
  n <- 5
  d <- 3
  pvalues <- c(0.4985223,
               0.1910987,
               0.5213573,
               0.1910409,
               0.6327676,
               0.6929403,
               0.6225058,
               0.8446667,
               0.8911884,
               0.2827703)
  set.seed(6)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   hybrid = FALSE,
                   cores = 1,
                   phenotypeCovariance = 'identity',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-7)
})


test_that("hybrid = FALSE. phenotypeCovariance = homogeneous", {
  # given
  p <- 10
  n <- 5
  d <- 3
  pvalues <- c(0.3591053,
               0.3617570,
               0.3562994,
               0.1201573,
               0.3377868,
               0.8581513,
               0.4901897,
               0.4350181,
               0.9916274,
               0.5053453)
  set.seed(6)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   hybrid = FALSE,
                   cores = 1,
                   phenotypeCovariance = 'homogeneous',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-7)
})