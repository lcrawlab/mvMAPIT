test_that("Simulate multiple phenotypes returns apropriate object", {
  # given
  p <- 10
  n <- 5
  d <- 3
  X <- matrix(runif(p * n), ncol = p)
  
  # when
  data <- simulate_phenotypes(X, d = d, logLevel = 'DEBUG')
  
  # then
  expect_equal(nrow(data$phenotype), n)
  expect_equal(ncol(data$phenotype), d)
  expect_equal(length(data), 7)
})
