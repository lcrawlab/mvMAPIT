/* Copyright 2017-2021 Lorin Crawford.
 *
 * This file uses the Catch unit testing library, alongside
 * testthat's simple bindings, to test a C++ function.
 *
 * For your own packages, ensure that your test files are
 * placed within the `src/` folder, and that you include
 * `LinkingTo: testthat` within your DESCRIPTION file.
 */
// probably need to make this an equivalent header file...
#include <RcppArmadillo.h>
#include "mapit/util.h"
#include "gsm/gsm.h"


// All test files should include the <testthat.h>
// header file.
#include <testthat.h>

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.
context("GetLinearKernel") {
  test_that("GetLinearKernel of identity is I/p") {
      arma::mat GSM = GetLinearKernel(arma::eye(2, 2));
    expect_true(arma::approx_equal(GSM, arma::eye(2, 2)/2, "absdiff", 0.001));
  }
}

context("ComputePCs") {
    // Problem from: https://mysite.science.uottawa.ca/phofstra/MAT2342/SVDproblems.pdf
    // However it has an error, U has wrong signs
    test_that("ComputePCs returns right values") {
        arma::mat X = { { 0, 1, 1 }, {sqrt(2), 2, 0}, {0, 1, 1} };
        arma::mat pc_cols = ComputePCs(X, 3);
        expect_true(pc_cols.n_cols == 3);
        arma::mat answer = { { -2/sqrt(3), sqrt(2)/sqrt(3), 0},
        {-4/sqrt(3), -sqrt(2)/sqrt(3), 0},
            {-2/sqrt(3), sqrt(2)/sqrt(3), 0} };
        expect_true(arma::approx_equal(pc_cols, answer, "absdiff", 0.01));
    }
}
