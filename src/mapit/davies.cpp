// Copyright 2021 Lorin Crawford.
#include "mapit/davies.h"
#include "mqs/mqs.h"

// #define WITH_LOGGER 1 // uncomment for logging during development
#ifdef WITH_LOGGER
#define SPDLOG_DISABLE_DEFAULT_LOGGER 1
#include <RcppSpdlog>
#endif

arma::mat davies_routine(const arma::mat &S, const arma::mat &Sinv,
                         const arma::mat &q,
                         const std::vector<arma::mat> &matrices) {
  arma::mat eigenvalues(matrices[0].n_cols, q.n_cols);

  for (int i = 0; i < q.n_cols; i++) {
    eigenvalues.col(i) = davies_routine_vec(S, Sinv, q.col(i), matrices);
  }
  return eigenvalues;
}

arma::vec davies_routine_vec(const arma::mat &S, const arma::mat &Sinv,
                             const arma::vec &q,
                             const std::vector<arma::mat> &matrices) {
  int num_variance_components = matrices.size();
  arma::vec evals;
  arma::vec eigval;
  arma::mat eigvec;

  arma::vec q_sub = q.subvec(1, num_variance_components - 1);
  arma::mat S_sub =
      S.submat(1, 1, num_variance_components - 1, num_variance_components - 1);
  arma::vec delta_null = arma::inv(S_sub) * q_sub;

  arma::mat A = arma::zeros(arma::size(matrices[0]));
  for (int i = 0; i < q_sub.n_elem; i++) {
    A = A + delta_null(i) * matrices[i + 1];
  }
  arma::eig_sym(eigval, eigvec, A);
  arma::mat EV_mat = compute_positive_ev_matrix(eigvec, eigval);
  evals = arma::eig_sym(EV_mat * compute_h_matrix(Sinv, matrices) * EV_mat);
  return evals;
}

arma::mat compute_positive_ev_matrix(const arma::mat &eigvec,
                                     const arma::mat &eigval) {
  arma::uvec ind_gt_zero = arma::find(eigval > 0);
  return eigvec.cols(ind_gt_zero) * arma::diagmat(sqrt(eigval(ind_gt_zero))) *
         arma::trans(eigvec.cols(ind_gt_zero));
}
