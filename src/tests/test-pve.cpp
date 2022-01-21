/* Copyright 2021 Lorin Crawford.
 */
#include "mapit/pve.h"
#include <RcppArmadillo.h>
#include <testthat.h>

context("compute_pve") {
  test_that("compute_pve for 3x3 variance components") {
    // given
    int component_index = 0;
    arma::mat A = arma::ones(4, 3);
    A(0, 1) = 2;
    A(0, 2) = 3;
    arma::rowvec correct_answer = {1.0 / 4, 2.0 / 5, 3.0 / 6};
    arma::rowvec correct_answer2 = {1.0 / 4, 1.0 / 5, 1.0 / 6};
    // when
    arma::rowvec result = compute_pve(A, component_index);
    arma::rowvec result2 = compute_pve(A, component_index + 1);
    // then
    expect_true(arma::approx_equal(result, correct_answer, "absdiff", 0.001));
    expect_true(arma::approx_equal(result2, correct_answer2, "absdiff", 0.001));
  }
}
