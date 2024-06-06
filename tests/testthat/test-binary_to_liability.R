test_that("binary_to_liability api test", {
  # given
  n_samples <- 10
  case_control_trait <- sample(0:1, n_samples, replace = TRUE)
  prevalence <- 0.1
  # when
  liabilities <- binary_to_liability(case_control_trait, prevalence)
  # then
  expect_true(is.numeric(liabilities))
  expect_equal(length(case_control_trait),
               length(liabilities))
})

test_that("binary_to_liability NA in case_control_trait", {
  # given
  n_samples <- 10
  n_missing <- 2
  case_control_trait <- sample(0:1, n_samples, replace = TRUE)
  prevalence <- 0.1
  case_control_trait[sample(1:n_samples, n_missing)] <- NA
  # when
  liabilities <- binary_to_liability(case_control_trait, prevalence)
  # then
  expect_equal(case_control_trait[is.na(case_control_trait)], liabilities[is.na(case_control_trait)])
  expect_equal(sum(!is.na(liabilities)), n_samples - n_missing)
  })