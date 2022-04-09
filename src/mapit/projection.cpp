// Copyright 2021 Lorin Crawford.
#include "mapit/projection.h"

arma::mat compute_projection_matrix(int n, const arma::mat &b) {
  arma::mat identity = arma::eye<arma::mat>(n, n);
  return identity - b * arma::inv(b.t() * b) * b.t();
}

arma::mat project_matrix(const arma::mat &A, const arma::mat &b) {
  arma::mat btb_inv = arma::inv(b.t() * b);
  return A - b * btb_inv * (b.t() * A) - (A * b) * btb_inv * b.t() +
         b * btb_inv * (b.t() * (A * b)) * btb_inv * b.t();
}
