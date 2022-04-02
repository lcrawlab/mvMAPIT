test_that("combinatorial test with d = 3. combinatorial.", {
  # given
  n <- 5
  d <- 3
  set.seed(853)
  phenotypes <- sprintf("p%s", 1:d)
  y <- matrix(runif(d * n), ncol = d)
  colnames(y) <- phenotypes
  correct_colnames <- c("p1*p1", "p2*p1", "p2*p2", "p3*p1", "p3*p2", "p3*p3")
  # when
  result <- mapit_struct_names(t(y))
  # then
  expect_equal(result, correct_colnames)
})

test_that("combinatorial test with d = 1. combinatorial.", {
  # given
  n <- 5
  d <- 1
  set.seed(853)
  phenotypes <- sprintf("p%s", 1:d)
  y <- matrix(runif(d * n), ncol = d)
  colnames(y) <- phenotypes
  correct_colnames <- c("p1")
  # when
  result <- mapit_struct_names(t(y))
  # then
  expect_equal(result, correct_colnames)
})

