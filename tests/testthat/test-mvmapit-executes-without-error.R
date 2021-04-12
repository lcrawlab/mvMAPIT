test_that("MvMapit executes without error when hybrid = TRUE.", {
  # given
  p <- 10
  n <- 5
  d <- 3
  X <- matrix(runif(p * n), n, p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X), t(Y), cores = 1, variantIndex = c(1))
  # then
  expect_equal(length(mapit$pvalues), p)
})

test_that("MvMapit executes without error when hybrid = FALSE.", {
  # given
  p <- 10
  n <- 5
  d <- 3
  X <- matrix(runif(p * n), n, p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X), t(Y), hybrid = FALSE, cores = 1, variantIndex = c(1))
  # then
  expect_equal(length(mapit$pvalues), p)
})

test_that("MvMapit executes without error when hybrid = FALSE and test = davies.", {
  # given
  p <- 10
  n <- 5
  d <- 3
  X <- matrix(runif(p * n), n, p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X), t(Y), hybrid = FALSE, test = 'davies', cores = 1, variantIndex = c(1))
  # then
  expect_equal(length(mapit$pvalues), p)
})
