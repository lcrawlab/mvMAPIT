// Copyright 2021 Lorin Crawford.
#include "mapit/projection.h"

arma::mat compute_projection_matrix(int n, const arma::mat &b) {
  arma::mat identity = arma::eye<arma::mat>(n, n);
  return identity - b * arma::inv(b.t() * b) * b.t();
}
