test_that("MvMapit executes without error when test = hybrid.", {
  # given
  p <- 10
  n <- 5
  d <- 1
  pvalues <- matrix(c(0.6772546,
               0.5690567,
               0.3366181,
               0.5523170,
               0.6307954,
               0.4436580,
               0.7090439,
               0.6147064,
               0.7128274,
               0.6088642), ncol = 1)
  set.seed(5)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   test = 'hybrid',
                   cores = 1,
                   accuracy = 1e-2,
                   phenotypeCovariance = 'covariance',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-7)
})

test_that("MvMapit executes without error when test = normal.", {
  # given
  p <- 10
  n <- 5
  d <- 3
  pvalues <- matrix(
              c(0.3409648,
               0.5879707,
               0.5973206,
               0.2912437,
               0.9626492,
               0.9496355,
               0.5612714,
               0.9947883,
               0.6359662,
               0.8031088), ncol = 1)
  set.seed(6)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   test = 'normal',
                   cores = 1,
                   phenotypeCovariance = 'covariance',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-7)
})

test_that("MvMapit executes without error when test = davies.", {
  # given
  p <- 10
  n <- 5
  d <- 3
  pvalues <- matrix(
             c(0.4855840,
               0.5042984,
               0.6235571,
               0.2560571,
               0.9150943,
               0.8862566,
               0.6406932,
               0.9294612,
               0.5974704,
               0.9704182), ncol = 1)
  set.seed(6)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   test = 'davies',
                   accuracy = 1e-5,
                   cores = 1,
                   variantIndex = c(1:p),
                   phenotypeCovariance = 'covariance',
                   logLevel = "ERROR")
  # then
  expect_equal(length(mapit$pvalues), p)
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-7)
})

test_that("MvMapit times computations when test = hybrid.", {
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
                   test = 'hybrid',
                   accuracy = 1e-2,
                   phenotypeCovariance = 'covariance',
                   logLevel = "ERROR")
  # then
  expect_equal(length(mapit$timings), 6)
})

test_that("MvMapit executes without error when C is not NULL, test = hybrid.", {
  # given
  p <- 10
  n <- 5
  d <- 1
  pvalues <- matrix(
             c(0.6828469,
               0.6014190,
               0.7815753,
               0.8516215,
               0.5315767,
               0.9810549,
               0.4784025,
               0.3926182,
               0.7893141,
               0.9439612), ncol = 1)
  set.seed(6)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  C <- matrix(runif(n * n), ncol = n)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   C = C,
                   cores = 1,
                   test = 'hybrid',
                   accuracy = 1e-6,
                   phenotypeCovariance = 'covariance',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-7)
})


test_that("test = normal. phenotypeCovariance = identity", {
  # given
  p <- 10
  n <- 5
  d <- 3
  pvalues <- matrix(
             c(0.4985223,
               0.1910987,
               0.5213573,
               0.1910409,
               0.6327676,
               0.6929403,
               0.6225058,
               0.8446667,
               0.8911884,
               0.2827703), ncol = 1)
  set.seed(6)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   cores = 1,
                   phenotypeCovariance = 'identity',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-7)
})

test_that("test = normal. phenotypeCovariance = homogeneous", {
  # given
  p <- 10
  n <- 5
  d <- 3
  pvalues <- matrix(
             c(0.3591053,
               0.3617570,
               0.3562994,
               0.1201573,
               0.3377868,
               0.8581513,
               0.4901897,
               0.4350181,
               0.9916274,
               0.5053453), ncol = 1)
  set.seed(6)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   accuracy = 1e-5,
                   cores = 1,
                   phenotypeCovariance = 'homogeneous',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-5)
})
