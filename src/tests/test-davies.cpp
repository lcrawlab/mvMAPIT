/* Copyright 2021 Lorin Crawford.
 */
#include "mapit/davies.h"
#include <RcppArmadillo.h>
#include <testthat.h>

context("davies_routine") {
  test_that("davies_routine for q matrix") {
    // given
    int num_samples = 4;
    int num_combinations = 6;
    int num_variance_components = 3;
    arma::mat q(num_variance_components, num_combinations);
    q.ones();
    std::vector<arma::mat> matrices(num_variance_components);
    arma::mat M(num_samples, num_samples);
    M.eye();
    arma::mat S(num_variance_components, num_variance_components);
    S.eye();
    matrices = {M, M, M};
    arma::mat correct_answer =
        arma::eye(num_combinations, num_variance_components);
    // when
    arma::mat result = davies_routine(S, S, q, matrices);
    // then
    expect_true(result.n_cols == num_combinations);
    expect_true(result.n_rows == num_samples);
  }
}
