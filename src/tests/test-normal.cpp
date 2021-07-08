/* Copyright 2021 Lorin Crawford.
 */
#include <RcppArmadillo.h>
#include <testthat.h>
#include "mapit/normal.h"

context("normal_pvalues") {
    test_that("normal_pvalues for 2x2") {
        // given
        int component_index = 0;
        arma::mat variance_estimate = arma::ones(2, 2);
        variance_estimate(0, 0) = 1.04;
        variance_estimate(1, 0) = 1.96;
        variance_estimate(0, 1) = 2.33;
        variance_estimate(1, 1) = 2.58;
        arma::mat standard_error = arma::ones(2, 2);
        arma::mat correct_answer = arma::ones(2, 2);
        correct_answer(0, 0) = 1 - 0.70;
        correct_answer(1, 0) = 1 - 0.95;
        correct_answer(0, 1) = 1 - 0.98;
        correct_answer(1, 1) = 1 - 0.99;
        // when
        arma::mat result = normal_pvalues(variance_estimate,
                                                    standard_error);
        // then
        expect_true(arma::approx_equal(result, correct_answer, "absdiff", 0.01));
    }
}
