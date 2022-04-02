test_that("pvalues and pve shows NA for variants not in variantIndex", {
  # given
  p <- 3
  n <- 10
  d <- 3
  set.seed(853)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  variantIndex <- c(1,3)
  otherIndex <- 2
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   test = 'normal',
                   cores = 1,
                   variantIndex = variantIndex,
                   logLevel = "DEBUG")
  # then
  expect_true(all(is.na(mapit$pves[otherIndex, ])))
  expect_true(all(is.na(mapit$pvalues[otherIndex, ])))
})

test_that("pve shows NA for covariance interactions", {
  # given
  p <- 3
  n <- 10
  d <- 3
  set.seed(853)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   test = 'normal',
                   cores = 1,
                   logLevel = "DEBUG")
  # then
  expect_true(all(is.na(mapit$pves[, 2])))
  expect_true(all(is.na(mapit$pves[, 4])))
  expect_true(all(is.na(mapit$pves[, 5])))
})

