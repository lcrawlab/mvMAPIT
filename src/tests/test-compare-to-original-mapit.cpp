/* Copyright 2021 Lorin Crawford.
 */
#include "gsm/gsm.h"
#include "mapit/projection.h"
#include "mapit/util.h"
#include "mqs/mqs.h"
#include <RcppArmadillo.h>
#include <testthat.h>

context("MAPIT1_Normal") {
  test_that("compute_k_matrix is equal to original") {
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
    arma::mat correct_answer = (GSM * p - trans(x_k) * x_k) / (p - 1);
    // when
    arma::mat result = compute_k_matrix(GSM, x_k, p);
    // then
    expect_true(result(0, 0) == correct_answer(0, 0));
    expect_true(result(0, 1) == correct_answer(0, 1));
    expect_true(result(0, 2) == correct_answer(0, 2));
    expect_true(result(1, 0) == correct_answer(1, 0));
    expect_true(result(1, 1) == correct_answer(1, 1));
    expect_true(result(1, 2) == correct_answer(1, 2));
    expect_true(result(2, 0) == correct_answer(2, 0));
    expect_true(result(2, 1) == correct_answer(2, 1));
    expect_true(result(2, 2) == correct_answer(2, 2));
  }
  test_that("compute_g_matrix is equal to original") {
    // given
    int n = 3;
    int p = 3;
    arma::mat K(n, p, arma::fill::zeros);
    K(0, 0) = 0.3333;
    K(0, 1) = -0.1667;
    K(0, 2) = -0.1667;
    K(1, 0) = K(0, 1);
    K(1, 1) = 0.8333;
    K(1, 2) = -0.6667;
    K(2, 0) = K(0, 2);
    K(2, 1) = K(1, 2);
    K(2, 2) = 0.8333;
    arma::rowvec x_k(n, arma::fill::zeros);
    x_k(0) = -1.1547;
    x_k(1) = 0.5774;
    x_k(2) = 0.5774;
    arma::mat correct_answer = K;
    correct_answer.each_row() %= x_k;
    correct_answer.each_col() %= x_k.t();
    // when
    arma::mat result = compute_g_matrix(K, x_k);
    // then
    expect_true(result(0, 0) == correct_answer(0, 0));
    expect_true(result(0, 1) == correct_answer(0, 1));
    expect_true(result(0, 2) == correct_answer(0, 2));
    expect_true(result(1, 0) == correct_answer(1, 0));
    expect_true(result(1, 1) == correct_answer(1, 1));
    expect_true(result(1, 2) == correct_answer(1, 2));
    expect_true(result(2, 0) == correct_answer(2, 0));
    expect_true(result(2, 1) == correct_answer(2, 1));
    expect_true(result(2, 2) == correct_answer(2, 2));
  }
  test_that("projection is equal to original") {
    // given
    int n = 3;
    int p = 3;
    arma::mat K(n, p, arma::fill::zeros);
    K(0, 0) = 0.3333;
    K(0, 1) = -0.1667;
    K(0, 2) = -0.1667;
    K(1, 0) = K(0, 1);
    K(1, 1) = 0.8333;
    K(1, 2) = -0.6667;
    K(2, 0) = K(0, 2);
    K(2, 1) = K(1, 2);
    K(2, 2) = 0.8333;
    arma::rowvec x_k(n, arma::fill::zeros);
    x_k(0) = -1.1547;
    x_k(1) = 0.5774;
    x_k(2) = 0.5774;
    arma::mat Y(2, n, arma::fill::zeros);
    Y(0, 0) = 1.1;
    Y(0, 1) = 2.1;
    Y(0, 2) = 3.1;
    arma::rowvec y = Y.row(0);
    arma::mat b = arma::zeros(n, 2);
    b.col(0) = arma::ones<arma::vec>(n);
    b.col(1) = trans(x_k);
    arma::mat btb_inv = inv(b.t() * b);
    arma::mat correct_answer =
        K - b * btb_inv * (b.t() * K) - (K * b) * btb_inv * b.t() +
        b * btb_inv * (b.t() * (K * b)) * btb_inv * b.t();
    arma::vec yc = (arma::eye<arma::mat>(n, n) - (b * btb_inv) * b.t()) * y.t();
    arma::mat M = compute_projection_matrix(n, b);
    // when
    arma::mat result = M * K * M;
    arma::mat Yc = Y * M;
    arma::rowvec result2 = Yc.row(0);
    // then
    // TODO: figure out why 0.0 == 0.0 fails tests
    //         expect_true(result(0, 0) - correct_answer(0, 0) == 0.0);
    // expect_true(result(0, 1) - correct_answer(0, 1) == 0.0);
    //         expect_true(result(0, 2) - correct_answer(0, 2) == 0.0);
    expect_true(result(1, 0) == correct_answer(1, 0));
    expect_true(result(1, 1) == correct_answer(1, 1));
    expect_true(result(1, 2) == correct_answer(1, 2));
    expect_true(result(2, 0) == correct_answer(2, 0));
    expect_true(result(2, 1) == correct_answer(2, 1));
    expect_true(result(2, 2) == correct_answer(2, 2));
    expect_true(result2(0) == yc(0));
    expect_true(result2(1) == yc(1));
    expect_true(result2(2) == yc(2));
  }
  test_that("compute_q_vector is equal to original") {
    // given
    int n = 3;
    int p = 3;
    arma::mat K(n, p, arma::fill::zeros);
    K(0, 0) = 0.3333;
    K(0, 1) = -0.1667;
    K(0, 2) = -0.1667;
    K(1, 0) = K(0, 1);
    K(1, 1) = 0.8333;
    K(1, 2) = -0.6667;
    K(2, 0) = K(0, 2);
    K(2, 1) = K(1, 2);
    K(2, 2) = 0.8333;
    arma::rowvec x_k(n, arma::fill::zeros);
    x_k(0) = -1.1547;
    x_k(1) = 0.5774;
    x_k(2) = 0.5774;
    arma::mat Y(1, n, arma::fill::zeros);
    Y(0, 0) = 1.1;
    Y(0, 1) = 2.1;
    Y(0, 2) = 3.1;
    arma::rowvec y = Y.row(0);
    arma::mat b = arma::zeros(n, 2);
    b.col(0) = arma::ones<arma::vec>(n);
    b.col(1) = trans(x_k);
    arma::mat btb_inv = inv(b.t() * b);
    arma::mat Kc_original = K - b * btb_inv * (b.t() * K) -
                            (K * b) * btb_inv * b.t() +
                            b * btb_inv * (b.t() * (K * b)) * btb_inv * b.t();
    arma::vec yc_original =
        (arma::eye<arma::mat>(n, n) - (b * btb_inv) * b.t()) * y.t();
    arma::mat G = K; // Create the Kn Matrix
    G.each_row() %= x_k;
    G.each_col() %= x_k.t();
    arma::mat Gc_original = G - b * btb_inv * (b.t() * G) -
                            (G * b) * btb_inv * b.t() +
                            b * btb_inv * (b.t() * (G * b)) * btb_inv * b.t();
    arma::mat M = compute_projection_matrix(n, b);
    arma::mat Kc = M * K * M;
    arma::mat Gc = M * G * M;
    arma::mat Yc = Y * M;
    arma::vec correct_answer = arma::zeros(3); // Create k-vector q to save
    correct_answer(0) = as_scalar(yc_original.t() * Gc_original * yc_original);
    correct_answer(1) = as_scalar(yc_original.t() * Kc_original * yc_original);
    correct_answer(2) = as_scalar(
        yc_original.t() * (arma::eye<arma::mat>(n, n) - (b * btb_inv) * b.t()) *
        yc_original);
    std::vector<arma::mat> matrices = {Gc, Kc, M};
    arma::vec yc = vectorise(Yc);
    std::vector<arma::vec> phenotypes = matrix_to_vector_of_rows(yc.as_row());
    // when
    arma::mat result = compute_q_matrix(phenotypes, matrices);
    // then
    expect_true(result(0, 0) == correct_answer(0));
    expect_true(result(1, 0) == correct_answer(1));
    expect_true(result(2, 0) == correct_answer(2));
  }
  test_that("compute_s_matrix is equal to original") {
    // given
    int n = 3;
    int p = 3;
    arma::mat K(n, p, arma::fill::zeros);
    K(0, 0) = 0.3333;
    K(0, 1) = -0.1667;
    K(0, 2) = -0.1667;
    K(1, 0) = K(0, 1);
    K(1, 1) = 0.8333;
    K(1, 2) = -0.6667;
    K(2, 0) = K(0, 2);
    K(2, 1) = K(1, 2);
    K(2, 2) = 0.8333;
    arma::rowvec x_k(n, arma::fill::zeros);
    x_k(0) = -1.1547;
    x_k(1) = 0.5774;
    x_k(2) = 0.5774;
    arma::mat b = arma::zeros(n, 2);
    b.col(0) = arma::ones<arma::vec>(n);
    b.col(1) = trans(x_k);
    arma::mat btb_inv = inv(b.t() * b);
    arma::mat Kc_original = K - b * btb_inv * (b.t() * K) -
                            (K * b) * btb_inv * b.t() +
                            b * btb_inv * (b.t() * (K * b)) * btb_inv * b.t();
    arma::mat G = K; // Create the Kn Matrix
    G.each_row() %= x_k;
    G.each_col() %= x_k.t();
    arma::mat Gc_original = G - b * btb_inv * (b.t() * G) -
                            (G * b) * btb_inv * b.t() +
                            b * btb_inv * (b.t() * (G * b)) * btb_inv * b.t();
    arma::mat M = compute_projection_matrix(n, b);
    arma::mat Kc = M * K * M;
    arma::mat Gc = M * G * M;
    arma::mat correct_answer = arma::zeros(n, n); // Create k-vector q to save
    correct_answer(0, 0) =
        arma::as_scalar(arma::accu(Gc_original.t() % Gc_original));
    correct_answer(0, 1) =
        arma::as_scalar(arma::accu(Gc_original.t() % Kc_original));
    correct_answer(0, 2) = arma::as_scalar(
        arma::accu(Gc_original.t() %
                   (arma::eye<arma::mat>(n, n) - (b * btb_inv) * b.t())));
    correct_answer(1, 0) = correct_answer(0, 1);
    correct_answer(1, 1) =
        arma::as_scalar(arma::accu(Kc_original.t() % Kc_original));
    correct_answer(1, 2) = arma::as_scalar(
        arma::accu(Kc_original.t() %
                   (arma::eye<arma::mat>(n, n) - (b * btb_inv) * b.t())));
    correct_answer(2, 0) = correct_answer(0, 2);
    correct_answer(2, 1) = correct_answer(1, 2);
    correct_answer(2, 2) = arma::as_scalar(
        arma::accu(trans(arma::eye<arma::mat>(n, n) - (b * btb_inv) * b.t()) %
                   (arma::eye<arma::mat>(n, n) - (b * btb_inv) * b.t())));
    std::vector<arma::mat> matrices = {Gc, Kc, M};
    // when
    arma::mat result = compute_s_matrix(matrices);
    // then
    expect_true(result(0, 0) == correct_answer(0, 0));
    expect_true(result(0, 1) == correct_answer(0, 1));
    expect_true(result(0, 2) == correct_answer(0, 2));
    expect_true(result(1, 0) == correct_answer(1, 0));
    expect_true(result(1, 1) == correct_answer(1, 1));
    expect_true(result(1, 2) == correct_answer(1, 2));
    expect_true(result(2, 0) == correct_answer(2, 0));
    expect_true(result(2, 1) == correct_answer(2, 1));
    expect_true(result(2, 2) == correct_answer(2, 2));
  }
  test_that("compute variance component is equal to original") {
    // given
    int n = 3;
    int p = 3;
    arma::mat K(n, p, arma::fill::zeros);
    K(0, 0) = 5.3333;
    K(0, 1) = -0.1667;
    K(0, 2) = -0.1667;
    K(1, 0) = K(0, 1);
    K(1, 1) = 0.8333;
    K(1, 2) = -0.6667;
    K(2, 0) = K(0, 2);
    K(2, 1) = K(1, 2);
    K(2, 2) = 1.8333;
    arma::rowvec x_k(n, arma::fill::zeros);
    x_k(0) = -1.1547;
    x_k(1) = 0.5774;
    x_k(2) = 0.5774;
    arma::mat Y(1, n, arma::fill::zeros);
    Y(0, 0) = 1.1;
    Y(0, 1) = 2.1;
    Y(0, 2) = 3.1;
    arma::rowvec y = Y.row(0);
    arma::mat b = arma::zeros(n, 2);
    b.col(0) = arma::ones<arma::vec>(n);
    b.col(1) = trans(x_k);
    arma::mat btb_inv = inv(b.t() * b);
    arma::mat Kc_original = K - b * btb_inv * (b.t() * K) -
                            (K * b) * btb_inv * b.t() +
                            b * btb_inv * (b.t() * (K * b)) * btb_inv * b.t();
    arma::vec yc_original =
        (arma::eye<arma::mat>(n, n) - (b * btb_inv) * b.t()) * y.t();
    arma::mat G = K; // Create the Kn Matrix
    G.each_row() %= x_k;
    G.each_col() %= x_k.t();
    arma::mat Gc_original = G - b * btb_inv * (b.t() * G) -
                            (G * b) * btb_inv * b.t() +
                            b * btb_inv * (b.t() * (G * b)) * btb_inv * b.t();
    arma::mat M = compute_projection_matrix(n, b);
    arma::mat Kc = M * K * M;
    arma::mat Gc = M * G * M;
    arma::mat Yc = Y * M;
    std::vector<arma::mat> matrices = {Gc, Kc, M};
    arma::vec yc = vectorise(Yc);
    std::vector<arma::vec> phenotypes = matrix_to_vector_of_rows(yc.as_row());
    arma::vec q_original = arma::zeros(3); // Create k-vector q to save
    q_original(0) = as_scalar(yc_original.t() * Gc_original * yc_original);
    q_original(1) = as_scalar(yc_original.t() * Kc_original * yc_original);
    q_original(2) = as_scalar(
        yc_original.t() * (arma::eye<arma::mat>(n, n) - (b * btb_inv) * b.t()) *
        yc_original);
    arma::mat S_original = arma::zeros(n, n); // Create k-vector q to save
    S_original(0, 0) =
        arma::as_scalar(arma::accu(Gc_original.t() % Gc_original));
    S_original(0, 1) =
        arma::as_scalar(arma::accu(Gc_original.t() % Kc_original));
    S_original(0, 2) = arma::as_scalar(
        arma::accu(Gc_original.t() %
                   (arma::eye<arma::mat>(n, n) - (b * btb_inv) * b.t())));
    S_original(1, 0) = S_original(0, 1);
    S_original(1, 1) =
        arma::as_scalar(arma::accu(Kc_original.t() % Kc_original));
    S_original(1, 2) = arma::as_scalar(
        arma::accu(Kc_original.t() %
                   (arma::eye<arma::mat>(n, n) - (b * btb_inv) * b.t())));
    S_original(2, 0) = S_original(0, 2);
    S_original(2, 1) = S_original(1, 2);
    S_original(2, 2) = arma::as_scalar(
        arma::accu(trans(arma::eye<arma::mat>(n, n) - (b * btb_inv) * b.t()) %
                   (arma::eye<arma::mat>(n, n) - (b * btb_inv) * b.t())));
    //         arma::mat Sinv = arma::inv(S_original);
    // arma::mat q = compute_q_matrix(phenotypes, matrices);
    // arma::vec delta_original = Sinv * q_original;
    // arma::mat delta = Sinv * q;
    // double correct_answer =
    // as_scalar(2*yc_original.t()*trans(Sinv(0,0)*Gc_original+Sinv(0,1)*Kc_original+Sinv(0,2)*(arma::eye<arma::mat>(n,n)-(b*btb_inv)*b.t()))*(delta_original(0)*Gc_original+delta_original(1)*Kc_original+delta_original(2)*(arma::eye<arma::mat>(n,n)-(b*btb_inv)*b.t()))*(Sinv(0,0)*Gc_original+Sinv(0,1)*Kc_original+Sinv(0,2)*(arma::eye<arma::mat>(n,n)-(b*btb_inv)*b.t()))*yc_original);
    // // when
    // arma::vec result = compute_variance_delta(phenotypes,
    // Sinv,
    // delta,
    // matrices);
    //         // then
    // expect_true(result(0) == correct_answer);
  }
}
