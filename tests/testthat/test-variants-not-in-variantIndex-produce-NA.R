test_that("pvalues and pve shows NA for variants not in variantIndex", {
  # given
  p <- 3
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
  variantIndex <- c(1,3)
  otherIndex <- 2
  # when
  mapit <- MvMAPIT(t(X),
                   t(Y),
                   test = 'normal',
                   cores = 1,
                   variantIndex = variantIndex,
                   phenotypeCovariance = 'combinatorial',
                   logLevel = "DEBUG")
  pves <- mapit$pves
  pves[!(c(1:nrow(pves)) %in% variantIndex)] <- NA
  # then
  expect_true(is.na(mapit$pves[otherIndex, 1]))
  expect_true(is.na(mapit$pvalues[otherIndex, 3]))
})
