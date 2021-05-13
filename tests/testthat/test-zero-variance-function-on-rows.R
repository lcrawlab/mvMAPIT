test_that("Zero variance does it's thing.", {
  # given
  p <- 10
  z <- 3
  o <- 2 * p - z
  zeros <- c(rep(0, z))
  ones <- c(rep(1, o))
  X <- matrix(c(zeros, ones), nrow = p)
  # when
  result <- remove_zero_variance(X)
  # then
  expect_equal(nrow(result), z)
})
