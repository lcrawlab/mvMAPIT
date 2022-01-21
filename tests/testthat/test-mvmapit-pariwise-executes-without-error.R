test_that("test = 'normal'. phenotypeCovariance = 'combinatorial'", {
  # given
  p <- 2
  n <- 10
  d <- 3
  pvalues <- matrix(
                    c(0.4705429,
                      0.5024536,
                      0.371318,
                      0.363676,
                      0.6591476,
                      0.6629176,
                      0.4778697,
                      0.6086180,
                      0.411222,
                      0.374,
                      0.0308,
                      0.05034315),
                    nrow = p, ncol = 6)
  colnames(pvalues) <- c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3")
  set.seed(853)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   test = 'normal',
                   cores = 1,
                   phenotypeCovariance = 'combinatorial',
                   logLevel = "DEBUG")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-4)
})

test_that("test = davies. phenotypeCovariance = 'combinatorial'", {
  # given
  p <- 2
  p <- 2
  n <- 10
  d <- 3
  pvalues <- matrix(
                    c(0.01558980,
                      0.01606852,
                      0.0,
                      0.0,
                      0.3795860,
                      0.3622172,
                      0.0002072584,
                      0.0962535911,
                      5.302820e-07,
                      2.003761e-05,
                      0.3353280,
                      0.2971691),
                    nrow = p, ncol = 6)
  colnames(pvalues) <- c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3")
  set.seed(853)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   test = 'davies',
                   cores = 1,
                   phenotypeCovariance = 'combinatorial',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-4)
})

test_that("phenotypeCovariance = 'combinatorial'", {
  # given
  p <- 2
  n <- 10
  d <- 3
  pvalues <- matrix(
                    c(0.4705429,
                      0.5024536,
                      0.371318,
                      0.363676,
                      0.6591476,
                      0.6629176,
                      0.4778697,
                      0.6086180,
                      0.411222,
                      0.374,
                      0.33532798,
                      0.05034315),
                    nrow = p, ncol = 6)
  colnames(pvalues) <- c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3")
  set.seed(853)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   test = 'hybrid',
                   cores = 1,
                   phenotypeCovariance = 'combinatorial',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-4)
})

test_that("C is not NULL. phenotypeCovariance = 'combinatorial'", {
  # given
  p <- 4
  n <- 10
  d <- 3
  pvalues <- matrix(c(0.7671534, 0.4078632, 0.6324734, 0.6599542, 0.5299002, 0.6711439,
                      0.6084251, 0.3368879, 0.2862572, 0.6728102, 0.5491810, 0.7753945,
                      0.4569704, 0.2566037, 0.4619346, 0.3094015, 0.7575085, 0.4011233,
                      0.3357123, 0.4324042, 0.7655489, 0.3377071, 0.8869128, 0.4731436),
                    nrow = p, ncol = 6, byrow = TRUE)
  colnames(pvalues) <- c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3")
  set.seed(199)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  C <- matrix(runif(n * n), ncol = n)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   C = C,
                   test = 'hybrid',
                   accuracy = 1e-5,
                   cores = 1,
                   phenotypeCovariance = 'combinatorial',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-4)
})

test_that("test = 'normal', C is not NULL. phenotypeCovariance = 'combinatorial'", {
  # given
  p <- 4
  n <- 10
  d <- 3
  pvalues <- matrix(c(0.7671534, 0.4078632, 0.6324734, 0.6599542, 0.5299002, 0.6711439,
                      0.6084251, 0.3368879, 0.2862572, 0.6728102, 0.5491810, 0.7753945,
                      0.4569704, 0.2566037, 0.4619346, 0.3094015, 0.7575085, 0.4011233,
                      0.3357123, 0.4324042, 0.7655489, 0.3377071, 0.8869128, 0.4731436),
                    nrow = p, ncol = 6, byrow = TRUE)
  colnames(pvalues) <- c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3")
  set.seed(199)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  C <- matrix(runif(n * n), ncol = n)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   C = C,
                   test = 'normal',
                   cores = 1,
                   phenotypeCovariance = 'combinatorial',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-4)
})

test_that("C is not NULL, test = 'davies'. phenotypeCovariance = 'combinatorial'", {
  # given
  p <- 4
  n <- 10
  d <- 3
  pvalues <- matrix(c(0.8475492, 0.0004668929, 0.97406290, 0.8097857, 0.6487072, 0.0000000,
                      0.7706045, 0.0041217079, 0.35731039, 0.6969462, 0.3551789, 0.1569766,
                      0.1958846, 0.8173460267, 0.08681906, 0.8106732, 0.2437403, 0.0000000,
                      0.0729410, 0.0000000000, 0.00000000, 0.7929982, 0.4209060, 0.7197998),
                    nrow = p, ncol = 6, byrow = TRUE)
  colnames(pvalues) <- c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3")
  set.seed(853)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  C <- matrix(runif(n * n), ncol = n)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   C = C,
                   test = 'davies',
                   accuracy = 1e-5,
                   cores = 1,
                   phenotypeCovariance = 'combinatorial',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-4)
})

test_that("test = 'davies'. phenotypeCovariance = 'combinatorial', d = 1", {
  # given
  p <- 10
  n <- 4
  d <- 1
  pvalues <- matrix(c(0.091008, 0.000524, 0.000000e+00, 0.840668, 0.000000e+00, 0.002915,
                      0.010521, 0.181760, 0.880265, 0.451950),
                    nrow = p, ncol = 1, byrow = TRUE)
  set.seed(853)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  C <- matrix(runif(n * n), ncol = n)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   C = C,
                   test = 'davies',
                   accuracy = 1e-5,
                   cores = 1,
                   phenotypeCovariance = 'combinatorial',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-4)
})
