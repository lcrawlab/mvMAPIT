test_that("MvMapit executes without error when hybrid = TRUE.", {
  # given
  p <- 10
  n <- 5
  d <- 1
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X), 
                   t(Y), 
                   cores = 1, 
                   #variantIndex = c(1),
                   phenotypeCovariance = 'covariance',
                   logLevel = "ERROR")
  # then
  expect_equal(length(mapit$pvalues), p)
})

test_that("MvMapit executes without error when hybrid = FALSE.", {
  # given
  p <- 10
  n <- 5
  d <- 3
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X), 
                   t(Y), 
                   hybrid = FALSE, 
                   cores = 1, 
                   variantIndex = c(1:2),
                   phenotypeCovariance = 'covariance', 
                   logLevel = "ERROR")
  # then
  expect_equal(length(mapit$pvalues), p)
})

test_that("MvMapit executes without error when hybrid = FALSE and test = davies.", {
  # given
  p <- 10
  n <- 5
  d <- 3
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X), 
                   t(Y), 
                   hybrid = FALSE, 
                   test = 'davies', 
                   cores = 1, 
                   variantIndex = c(1:2),
                   phenotypeCovariance = 'covariance',
                   logLevel = "ERROR")
  # then
  expect_equal(length(mapit$pvalues), p)
})

test_that("MvMapit times computations when hybrid = TRUE.", {
  # given
  p <- 10
  n <- 5
  d <- 1
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   cores = 1,
                   #variantIndex = c(1),
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
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  C <- matrix(runif(n * n), ncol = n)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   C = C,
                   cores = 1,
                   #variantIndex = c(1),
                   hybrid = TRUE,
                   phenotypeCovariance = 'covariance',
                   logLevel = "ERROR")
  # then
  expect_equal(length(mapit$pvalues), p)
})
