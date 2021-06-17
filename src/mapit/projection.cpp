// Copyright 2021 Lorin Crawford.
#include <RcppArmadillo.h>
#include "projection.h"

arma::mat ComputeProjectionMatrix(int n, arma::mat b) {
    arma::mat identity = arma::eye<arma::mat>(n, n);
    return identity - b * arma::inv(b.t() * b) * b.t();
}
