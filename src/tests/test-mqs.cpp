/* Copyright 2021 Lorin Crawford.
 */
#include <RcppArmadillo.h>
#include <testthat.h>
#include "mqs/mqs.h"

context("product_trace") {
    test_that("product_trace 2x2 identity") {
        // given
        arma::mat A = arma::eye(2, 2);
        arma::mat B = arma::eye(2, 2);
        double correct_answer = 2;
        // when
        double result = product_trace(A, B);
        // then
        expect_true(result == correct_answer);
    }
    test_that("product_trace 2x2 ones") {
        // given
        arma::mat A = arma::ones(2, 2);
        arma::mat B = arma::ones(2, 2);
        double correct_answer = 4;
        // when
        double result = product_trace(A, B);
        // then
        expect_true(result == correct_answer);
    }
}

context("compute_q_matrix") {
    test_that("compute_q_matrix for 3 phenotypes returns 6 combinations") {
        // given
        int num_samples = 3;
        int num_phenotypes = 3;
        int num_variance_components = 3;
        int num_combinations = 6;
        std::vector<arma::vec> V(num_phenotypes);
        std::vector<arma::mat> M(num_variance_components);
        arma::mat matrix(num_samples, num_samples); matrix.eye();
        arma::vec v1(num_samples, arma::fill::zeros); v1(0) = 1;
        arma::vec v2(num_samples, arma::fill::zeros); v2(1) = 1;
        arma::vec v3(num_samples, arma::fill::zeros); v3(2) = 1;
        V[0] = v1;
        V[1] = v2;
        V[2] = v3;
        M = { matrix, matrix, matrix };
        // when
        arma::mat result = compute_q_matrix(V, M);
        // then
        expect_true(result.n_rows == num_variance_components);
        expect_true(result.n_cols == num_combinations);
    }
    test_that("compute_q_matrix for 2 phenotypes returns 3 combinations") {
        // given
        int num_samples = 3;
        int num_phenotypes = 2;
        int num_variance_components = 3;
        int num_combinations = 3;
        std::vector<arma::vec> V(num_phenotypes);
        std::vector<arma::mat> M(num_variance_components);
        arma::mat matrix(num_samples, num_samples); matrix.eye();
        arma::vec v1(num_samples, arma::fill::zeros); v1(0) = 1;
        arma::vec v2(num_samples, arma::fill::zeros); v2(1) = 1;
        V[0] = v1;
        V[1] = v2;
        M = { matrix, matrix, matrix };
        // when
        arma::mat result = compute_q_matrix(V, M);
        // then
        expect_true(result.n_rows == num_variance_components);
        expect_true(result.n_cols == num_combinations);
    }
}

context("compute_q_vector") {
    test_that("compute_q_vector identity") {
        // given
        int num_samples = 3;
        int num_variance_components = 3;
        int num_combinations = 6;
        std::vector<arma::mat> M(num_variance_components);
        arma::mat matrix(num_samples, num_samples); matrix.eye();
        arma::vec v1(num_samples, arma::fill::ones);
        arma::vec v2(num_samples, arma::fill::ones);
        M = { matrix, 2 * matrix, 3 * matrix };
        arma::vec correct_answer(num_variance_components);
        correct_answer(0) = 3;
        correct_answer(1) = 6;
        correct_answer(2) = 9;
        // when
        arma::vec result = compute_q_vector(v1, v2, M);
        // then
        //expect_true(arma::approx_equal(result, correct_answer, "absdiff", 0.001));
        expect_true(result(0) == correct_answer(0));
        expect_true(result(1) == correct_answer(1));
        expect_true(result(2) == correct_answer(2));
    }
    test_that("compute_q_vector homogeneous") {
        // given
        int num_samples = 3;
        int num_variance_components = 3;
        int num_combinations = 6;
        std::vector<arma::mat> M(num_variance_components);
        arma::mat matrix(num_samples, num_samples); matrix.ones();
        arma::vec v1(num_samples, arma::fill::ones);
        arma::vec v2(num_samples, arma::fill::ones);
        M = { matrix, matrix, matrix };
        arma::vec correct_answer(num_variance_components);
        correct_answer(0) = 9;
        correct_answer(1) = 9;
        correct_answer(2) = 9;
        // when
        arma::vec result = compute_q_vector(v1, v2, M);
        // then
        //expect_true(arma::approx_equal(result, correct_answer, "absdiff", 0.001));
        expect_true(result(0) == correct_answer(0));
        expect_true(result(1) == correct_answer(1));
        expect_true(result(2) == correct_answer(2));
    }
}
