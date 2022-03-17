test_that("test = 'normal'. phenotypeCovariance = 'combinatorial'", {
  # given
  p <- 2
  n <- 10
  d <- 3
  pvalues <- matrix(c(
0.4990573, 0.4478648, 0.9574136, 0.4662016, 0.4782672, 0.1381317, 0.612,
0.5015375, 0.4619467, 0.1347061, 0.4507170, 0.6405290, 0.2410251, 0.425),
                    nrow = p, ncol = 7, byrow = TRUE)
  colnames(pvalues) <- c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3", "metap")
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
  n <- 10
  d <- 3
  pvalues <- matrix(c(
0.01624319,     0, 0.6531582, 1.380850e-04, 0.000000e+00, 0.419213, 0.0,
0.02694842,     0, 0.3977345, 5.755901e-05, 1.216344e-08, 0.490887, 0.0),
                    nrow = p, ncol = 7, byrow = TRUE)
  colnames(pvalues) <- c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3", "metap")
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
  pvalues <- matrix(c(
0.4990573, 0.4478648, 0.9574136, 0.4662016, 0.4782672, 0.1381317, 0.612,
0.5015375, 0.4619467, 0.1347061, 0.4507170, 0.6405290, 0.2410251, 0.425),
                    nrow = p, ncol = 7, byrow = TRUE)
  colnames(pvalues) <- c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3", "metap")
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
0.6876487, 0.2148062, 0.5640931, 0.1657485, 0.2837563, 0.5020969, 0.409,
0.8920097, 0.9107812, 0.9787608, 0.6248188, 0.2751130, 0.4994958, 0.945,
0.5868067, 0.5128342, 0.3823134, 0.8747280, 0.2352273, 0.6889640, 0.767,
0.3184337, 0.5047131, 0.6774045, 0.3071930, 0.8257162, 0.5527816, 0.756
                      ),
                    nrow = p, ncol = 7, byrow = TRUE)
  colnames(pvalues) <- c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3", "metap")
  set.seed(29)
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
  pvalues <- matrix(c(
0.6876487, 0.2148062, 0.5640931, 0.1657485, 0.2837563, 0.5020969, 0.409,
0.8920097, 0.9107812, 0.9787608, 0.6248188, 0.2751130, 0.4994958, 0.945,
0.5868067, 0.5128342, 0.3823134, 0.8747280, 0.2352273, 0.6889640, 0.767,
0.3184337, 0.5047131, 0.6774045, 0.3071930, 0.8257162, 0.5527816, 0.756
                      ),
                    nrow = p, ncol = 7, byrow = TRUE)
  colnames(pvalues) <- c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3", "metap")
  set.seed(29)
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
  pvalues <- matrix(c(
0.13395441,     0,     0, 0.8361945, 0.001629745, 0.009014612, 0,
0.00000000,     0,     0, 0.5526955, 0.000000000, 0.000000000, 0,
0.01039137,     0,     0, 0.7475025, 0.000000000, 0.000000000, 0,
0.00000000,     0,     0, 0.7418526, 0.000000000, 0.000000000, 0
),
                    nrow = p, ncol = 7, byrow = TRUE)
  colnames(pvalues) <- c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3", "metap")
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
  pvalues <- matrix(c(0.0000000,
                0.7422639,
                0.7635732,
                0.0000000,
                0.1059719,
                0.5984818,
                0.7906796,
                0.0000000,
                0.0000000,
                0.4871070),
                    nrow = p, ncol = 1, byrow = TRUE)
  set.seed(20)
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
