/* Copyright 2021 Lorin Crawford.
 */
#include <RcppArmadillo.h>
#include <testthat.h>
#include "mapit/pve.h"

context("compute_pve") {
    test_that("compute_pve for 3x3 variance components") {
        // given
        arma::mat A = arma::ones(3, 3);
        A(1, 0) = 2;
        A(2, 0) = 3;
        arma::vec correct_answer(3);
        correct_answer(0) = 1 / 2;
        correct_answer(1) = 2 / 4;
        correct_answer(2) = 3 / 5;
        // when
        arma::rowvec result = compute_pve(A);
        // then
        expect_true(result(0) == correct_answer(1) );
        expect_true(result(2) == correct_answer(0) );
        expect_true(result(1) == correct_answer(2) );
        expect_true(arma::approx_equal(result, correct_answer, "absdiff", 0.001));
    }
}
