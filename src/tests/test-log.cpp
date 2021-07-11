/* Copyright 2021 Lorin Crawford.
 */
#include <RcppArmadillo.h>
#include <testthat.h>
#include "logging/log.h"

context("vector_to_string") {
    test_that("vector_to_string for 3 elements") {
        // given
        arma::vec v = {1.0, 2.5, 3.0};
        std::string correct_answer = "1.000000, 2.500000, 3.000000";
        // when
        std::string result = vector_to_string(v);
        // then
        expect_true(result == correct_answer);
    }
    test_that("matrix_to_string for 3x3 elements") {
        // given
        arma::mat M(3, 3); M.ones();
        std::string correct_answer =
                    "1.000000, 1.000000, 1.000000, \n"
                    "1.000000, 1.000000, 1.000000, \n"
                    "1.000000, 1.000000, 1.000000, \n";
        // when
        std::string result = matrix_to_string(M);
        // then
        expect_true(result == correct_answer);
    }
}
