test_that("hybrid = FALSE. phenotypeCovariance = ''", {
  # given
  p <- 2
  n <- 10
  d <- 3
  pvalues <- matrix(c(0.4705, 0.5025, 0.4007, 0.3927, 0.6591, 0.6629, 0.4171, 0.0285, 0.3101, 0.1738, 0.0308, 0.0503),
                    nrow = p, ncol = 6)
  colnames(pvalues) <- c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3")
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
                      0.09625359,
                      0.3101424,
                      0.1737842,
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
                   hybrid = TRUE,
                   cores = 1,
                   phenotypeCovariance = '',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-4)
})

test_that("hybrid = TRUE, C is not NULL. phenotypeCovariance = ''", {
  # given
  p <- 4
  n <- 10
  d <- 3
  pvalues <- matrix(c(0.1659513, 0.2242929, 0.4392573, 0.4798231, 0.5942647, 0.9670782,
                      0.6098269, 0.2332077, 0.6874865, 0.2094195, 0.6129602, 0.6231966,
                      0.8408090, 0.1349153, 0.8036961, 0.2815634, 0.9890769, 0.3842330,
                      0.1122040, 0.4313260, 0.2153343, 0.8230605, 0.4501949, 0.8414133),
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
                   hybrid = TRUE,
                   accuracy = 1e-5,
                   cores = 1,
                   phenotypeCovariance = '',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-4)
})

test_that("hybrid = FALSE, C is not NULL. phenotypeCovariance = ''", {
  # given
  p <- 4
  n <- 10
  d <- 3
  pvalues <- matrix(c(0.1659513, 0.6098269, 0.8408090, 0.1122040, 0.2242929, 0.2332077,
                      0.1349153, 0.4313260, 0.4392573, 0.6874865, 0.8036961, 0.2153343,
                      0.4798231, 0.2094195, 0.2815634, 0.8230605, 0.5942647, 0.6129602,
                      0.9890769, 0.4501949, 0.9670782, 0.6231966, 0.3842330, 0.8414133),
                    nrow = p, ncol = 6)
  colnames(pvalues) <- c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3")
  set.seed(853)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  C <- matrix(runif(n * n), ncol = n)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   C = C,
                   hybrid = FALSE,
                   cores = 1,
                   phenotypeCovariance = '',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-4)
})

test_that("hybrid = FALSE, C is not NULL, test = 'davies'. phenotypeCovariance = ''", {
  # given
  p <- 4
  n <- 10
  d <- 3
  pvalues <- matrix(c(0.8920589, 0.0, 0.8704164, 0.9564335, 0.126598227, 0.92516224,
                      0.6506412, 5.062263e-07, 0.5693004, 0.5624713, 0.001549884, 0.53344421,
                      0.7103750, 1.108789e-04, 0.8424372, 0.6119419, 0.977816063, 0.04392313,
                      0.1437574, 0.5462071, 0.0, 0.8068312, 0.571373467, 0.87882231),
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
                   hybrid = FALSE,
                   test = 'davies',
                   accuracy = 1e-8,
                   cores = 1,
                   phenotypeCovariance = '',
                   logLevel = "ERROR")
  # then
  print(mapit$pvalues)
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-4)
})



test_that("hybrid = FALSE, test = 'davies'. phenotypeCovariance = '', d = 1", {
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
                   hybrid = FALSE,
                   test = 'davies',
                   accuracy = 1e-5,
                   cores = 1,
                   phenotypeCovariance = '',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-4)
})
