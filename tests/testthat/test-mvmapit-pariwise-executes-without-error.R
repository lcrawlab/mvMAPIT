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
                    c(0.01597664,
                      0.01640104,
                      0.0,
                      0.0,
                      0.3799469,
                      0.3618748,
                      0.0003248277,
                      0.0966434584,
                      0.0,
                      0.0,
                      0.3355311,
                      0.2968991),
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
                      0.09664346,
                      0.3101424,
                      0.1737842,
                      0.33553110,
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

test_that("hybrid = TRUE, C is not NULL. phenotypeCovariance = ''", {
  # given
  p <- 4
  n <- 10
  d <- 3
  pvalues <- matrix(c(0.1659513, 0.6098269, 0.8408090, 0.1122040, 0.2242929, 0.2332077,
                      0.1349153, 0.4313260, 0.4392573, 0.6874865, 0.8036961, 0.2153343,
                      0.4798231, 0.2094195, 0.2815634, 0.8230605, 0.5942647, 0.6129602,
                      0.9890769, 0.4501949, 0.9670782, 0.6231966, 0.3842330, 0.8414133),
                    nrow = p, ncol = 6)
  set.seed(853)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  C <- matrix(runif(n * n), ncol = n)
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   C = C,
                   hybrid = TRUE,
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
  pvalues <- matrix(c(0.8917604, 0.0, 0.8697775, 0.9576408, 0.126735221, 0.9264005,
                      0.6499192, 0.0, 0.5692051, 0.5625198, 0.001621203, 0.5336411,
                      0.7108441, 0.0001347479, 0.8417294, 0.6122342, 0.978414917, 0.0438949,
                      0.1439900, 0.5461373817, 0.0, 0.8058156, 0.571396524, 0.8780500),
                    nrow = p, ncol = 6, byrow = TRUE)
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
                   cores = 1,
                   phenotypeCovariance = '',
                   logLevel = "ERROR")
  # then
  expect_equal(mapit$pvalues, pvalues, tolerance = 1e-4)
})
