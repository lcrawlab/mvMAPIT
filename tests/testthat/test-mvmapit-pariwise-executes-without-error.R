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
  pvalues <- matrix(c(
0.177918, 0.115599, 0.620124, 0.150483, 0.598889, 0.7013684,
0.574333, 0.852303, 0.282658, 0.398936, 0.260208, 0.5625580,
0.703179, 0.332805, 0.345765, 0.292213, 0.634062, 0.5314551,
0.212619, 0.195149, 0.542377, 0.285747, 0.285560, 0.1112160),
                    nrow = p, ncol = 6, byrow = TRUE)
  colnames(pvalues) <- c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3")
  set.seed(23453)
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
  pvalues <- matrix(c(0.5837010, 0.6454003, 0.7744985, 0.5666703, 0.6525211, 0.5808023,
                      0.5283121, 0.9156183, 0.7902189, 0.5965195, 0.8322132, 0.2814996,
                      0.8475289, 0.6402617, 0.4896582, 0.9401121, 0.4925769, 0.5610111,
                      0.9140411, 0.4819994, 0.7677112, 0.3843962, 0.6217984, 0.4187421),
                    nrow = p, ncol = 6, byrow = TRUE)
  colnames(pvalues) <- c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3")
  set.seed(82853)
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
  pvalues <- matrix(c(0.8920589, 0.0, 0.8704164, 0.64820, 0.21871, 0.92516224,
                      0.6506412, 0.00109, 0.5693004, 0.33923, 0.00889, 0.53344421,
                      0.7103750, 0.81909, 0.8424372, 0.92670, 0.517, 0.04392313,
                      0.1437574, 0.80281, 0.0, 0.96080, 0.6337, 0.87882231),
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
  pvalues <- matrix(c(7.097669e-01,7.849280e-07, 0.000000e+00, 8.378058e-01, 0.000000e+00, 1.899326e-01,
                      4.319290e-01, 5.793326e-01, 9.784368e-01, 3.942235e-01),
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
