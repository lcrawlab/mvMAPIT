test_that("pairwise test with d = 3. test = 'normal'.", {
  # given
  p <- 10
  n <- 5
  d <- 3
  pvalues <- matrix(c(0.4705, 0.5025, 0.4007, 0.3927, 0.6591, 0.6629, 0.4171, 0.0285, 0.3101, 0.1738, 0.0308, 0.0503),
                    nrow = p, ncol = 6)
  set.seed(853)
  variants <- sprintf("SNP%s", 1:p)
  phenotypes <- sprintf("Q%s", 1:d)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  colnames(X) <- variants
  colnames(Y) <- phenotypes
  resulting_colnames <- c("Q1*Q1", "Q2*Q1", "Q2*Q2", "Q3*Q1", "Q3*Q2", "Q3*Q3")
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   test = 'normal',
                   cores = 1,
                   phenotypeCovariance = 'combinatorial',
                   logLevel = "ERROR")
  # then
  # print(mapit$pvalues)
  expect_equal(rownames(mapit$pvalues), variants)
  expect_equal(colnames(mapit$pvalues), resulting_colnames)
})

test_that("pairwise test with d = 3. test = 'normal'.", {
  # given
  p <- 10
  n <- 5
  d <- 3
  pvalues <- matrix(c(0.4705, 0.5025, 0.4007, 0.3927, 0.6591, 0.6629, 0.4171, 0.0285, 0.3101, 0.1738, 0.0308, 0.0503),
                    nrow = p, ncol = 6)
  set.seed(853)
  variants <- sprintf("SNP%s", 1:p)
  phenotypes <- sprintf("P%s", 1:d)
  X <- matrix(runif(p * n), ncol = p)
  Y <- matrix(runif(d * n), ncol = d)
  colnames(X) <- variants
  colnames(Y) <- phenotypes
  resulting_colnames <- c("kronecker")
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   test = 'normal',
                   cores = 1,
                   phenotypeCovariance = 'identity',
                   logLevel = "ERROR")
  # then
  # print(mapit$pvalues)
  expect_equal(rownames(mapit$pvalues), variants)
  expect_equal(colnames(mapit$pvalues), resulting_colnames)
})
