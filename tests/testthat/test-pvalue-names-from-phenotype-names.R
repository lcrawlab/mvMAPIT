test_that("combinatorial test with d = 3. combinatorial.", {
  # given
  phenotypeCovariance <- 'combinatorial'
  n <- 5
  d <- 3
  set.seed(853)
  phenotypes <- sprintf("P%s", 1:d)
  Y <- matrix(runif(d * n), ncol = d)
  colnames(Y) <- phenotypes
  correct_colnames <- c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3")
  # when
  result <- mapit_struct_names(t(Y), phenotypeCovariance)
  # then
  expect_equal(result, correct_colnames)
})

test_that("combinatorial test with d = 3. Kronecker.", {
  # given
  phenotypeCovariance <- 'identity'
  n <- 5
  d <- 3
  set.seed(853)
  phenotypes <- sprintf("P%s", 1:d)
  Y <- matrix(runif(d * n), ncol = d)
  colnames(Y) <- phenotypes
  correct_colnames <- c("kronecker")
  # when
  result <- mapit_struct_names(t(Y), phenotypeCovariance)
  # then
  expect_equal(result, correct_colnames)
})
