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
