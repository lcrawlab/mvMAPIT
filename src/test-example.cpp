/* Copyright 2017-2021 Lorin Crawford.
 *
 * This file uses the Catch unit testing library, alongside
 * testthat's simple bindings, to test a C++ function.
 *
 * For your own packages, ensure that your test files are
 * placed within the `src/` folder, and that you include
 * `LinkingTo: testthat` within your DESCRIPTION file.
 */
#include <RcppArmadillo.h>


// All test files should include the <testthat.h>
// header file.
#include <testthat.h>

#include "mapit/util.h"
#include "gsm/gsm.h"

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.
context("get_linear_kernel") {
  test_that("get_linear_kernel of identity is I/p") {
      // given
      arma::mat I = arma::eye(2, 2);
      // when
      arma::mat GSM = get_linear_kernel(I);
      // then
      expect_true(arma::approx_equal(GSM, arma::eye(2, 2)/2, "absdiff", 0.001));
  }
}

context("compute_principal_components") {
    // Problem from: https://mysite.science.uottawa.ca/phofstra/MAT2342/SVDproblems.pdf
    // However it has an error, U has wrong signs
    test_that("compute_principal_components returns right values") {
        // given
        arma::mat X = { { 0, 1, 1 }, {sqrt(2), 2, 0}, {0, 1, 1} };
        arma::mat answer = { { -2/sqrt(3), sqrt(2)/sqrt(3), 0},
        {-4/sqrt(3), -sqrt(2)/sqrt(3), 0},
            {-2/sqrt(3), sqrt(2)/sqrt(3), 0} };
        // when
        arma::mat pc_cols = compute_principal_components(X, 3);
        // then
        expect_true(pc_cols.n_cols == 3);
        expect_true(arma::approx_equal(pc_cols, answer, "absdiff", 0.01));
    }
}

context("skip_variant") {
    test_that("test_variant does not skip when empty") {
        // given
        int i = 1;
        arma::vec ind;
        // when
        bool skip = skip_variant(ind, i);
        // then
        expect_true(skip == false);
    }
    test_that("test_variant skips when not in list") {
        // given
        int i = 10;
        arma::vec ind(1);  ind = 2;
        // when
        bool skip = skip_variant(ind, i);
        // then
        expect_true(skip == true);
    }
    test_that("test_variant does not skip when in list") {
        // given
        int i = 1;
        arma::vec ind(1);  ind = i + 1;
        // when
        bool skip = skip_variant(ind, i);
        // then
        expect_true(skip == false);
    }
}

context("matrix_to_vector_of_vectors") {
    test_that("matrix_to_vector_of_vectors returns correct vectors") {
        // given
        arma::mat matrix(3, 3); matrix.eye();
        arma::vec v1(3, arma::fill::zeros); v1(0) = 1;
        arma::vec v2(3, arma::fill::zeros); v2(1) = 1;
        arma::vec v3(3, arma::fill::zeros); v3(2) = 1;
        // when
        std::vector<arma::vec> vectors = matrix_to_vector_of_vectors(matrix);
        // then
        expect_true(arma::approx_equal(vectors[0], v1, "absdiff", 0.01));
        expect_true(arma::approx_equal(vectors[1], v2, "absdiff", 0.01));
        expect_true(arma::approx_equal(vectors[2], v3, "absdiff", 0.01));
    }
}
