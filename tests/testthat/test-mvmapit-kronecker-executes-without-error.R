test_that("MvMapit executes without error when test = hybrid.", {
  # given
  p <- 10
  n <- 5
  d <- 1
  pvalues <- matrix(c(
0.7315409,
0.7961257,
0.7607653,
0.4380301,
0.5265725,
0.8722502,
0.5097280,
0.6620606,
0.4780134,
0.8188633
), ncol = 1)
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
  pvalues <- matrix(c(
0.5823842,
0.4968245,
0.9075542,
0.8032449,
0.6508556,
0.6345080,
0.2701859,
0.7694585,
0.8733472,
0.7751935
), ncol = 1)
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
  pvalues <- matrix(c(
0.5649081,
0.5902518,
0.9844545,
0.5968390,
0.3404326,
0.6872187,
0.3451383,
0.4676168,
0.9013817,
0.6097202
), ncol = 1)
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
  pvalues <- matrix(c(
0.6627869,
0.5757486,
0.4161929,
0.3638452,
0.9156424,
0.6960296,
0.7774770,
0.5700789,
0.7186994,
0.6525567
), ncol = 1, byrow = TRUE)
  set.seed(23)
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
  pvalues <- matrix(c(
0.4852139,
0.5064176,
0.8818390,
0.4643340,
0.5941437,
0.7743915,
0.2540507,
0.8820358,
0.3711353,
0.2102760
               ), ncol = 1)
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
  pvalues <- matrix(c(
0.8384163,
0.8715925,
0.7138586,
0.7740888,
0.5336523,
0.4166844,
0.8883106,
0.9463360,
0.4320256,
0.4812932), ncol = 1, byrow = TRUE)
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
