/* Copyright 2022 Julian Stamp.
 */
#include <chrono>
#include <iostream>
#include <vector>
#include "mapit/projection.h"
#include "gsm/gsm.h"
#include <RcppArmadillo.h>
#include <testthat.h>

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

context("projection as vector multiplication") {
  test_that("projection for small matrix") {
    // given
    int n = 3;
    int p = 3;
    arma::mat X(n, p, arma::fill::zeros);
    X(0, 0) = -1.1547;
    X(0, 1) = 0.5774;
    X(0, 2) = 0.5774;
    X(1, 0) = -0.5774;
    X(1, 1) = -0.5774;
    X(1, 2) = 1.1547;
    X(2, 0) = -0.5774;
    X(2, 1) = 1.1547;
    X(2, 2) = -0.5774;
    arma::rowvec x_k(n, arma::fill::zeros);
    x_k(0) = -1.1547;
    x_k(1) = 0.5774;
    x_k(2) = 0.5774;
    arma::mat GSM = get_linear_kernel(X);
    arma::mat b = arma::zeros(n, 2);
    b.col(0) = arma::ones<arma::vec>(n);
    b.col(1) = arma::trans(x_k);
    arma::mat M = compute_projection_matrix(n, b);
    arma::mat K = compute_k_matrix(GSM, x_k, p);
    arma::mat correct_answer = M * K * M;
    // when
    arma::mat result = project_matrix(K, b);
    // then
    expect_true(arma::approx_equal(result, correct_answer, "absdiff", 0.00001));
    // expect_true(result(0, 0) == correct_answer(0, 0)); // -0.0 == 0.0 error
    // expect_true(result(0, 1) == correct_answer(0, 1));
    // expect_true(result(0, 2) == correct_answer(0, 2));
    // expect_true(result(1, 0) == correct_answer(1, 0));
    expect_true(result(1, 1) == correct_answer(1, 1));
    expect_true(result(1, 2) == correct_answer(1, 2));
    // expect_true(result(2, 0) == correct_answer(2, 0));
    expect_true(result(2, 1) == correct_answer(2, 1));
    expect_true(result(2, 2) == correct_answer(2, 2));
  }
}

// Uncomment if want to run speed test.
// context("speed test projection") {
  // test_that("matrix multiplication vs project_matrix") {
    // // given
    // int l = 2;
    // int k = 100;
    // int N[] = {1000, 2000};
    // arma::mat execution_t(k, 2);
    // arma::mat mean_time(l, 2);
    // for (int i = 0; i < l; i++) {
      // std::cout << N[i] << "\n";
      // arma::mat X;
      // X.randu(N[i], N[i]);
      // arma::mat b = arma::zeros(N[i], 2);
      // b.col(0) = arma::ones<arma::vec>(N[i]);
      // b.col(1) = arma::vec(N[i], arma::fill::randu);
    // // when
      // for (int j = 0; j < k; j++) {
        // auto start = high_resolution_clock::now();
        // arma::mat R = X * X * X;
        // auto end = high_resolution_clock::now();
        // execution_t(j, 0) = duration_cast<milliseconds>(end - start).count();
        // start = high_resolution_clock::now();
        // R = project_matrix(X, b);
        // end = high_resolution_clock::now();
        // execution_t(j, 1) = duration_cast<milliseconds>(end - start).count();
      // }
      // mean_time.row(i) = arma::mean(execution_t, 0);
      // std::cout << mean_time.row(i) << "\n";
    // // then
      // expect_true(mean_time(i, 0) > mean_time(i, 1));
    // }
  // }
// }
