// Copyright 2021 Lorin Crawford.
#include <RcppArmadillo.h>
#include "gsm/gsm.h"

arma::mat get_linear_kernel(const arma::mat& X) {
    double p = X.n_rows;
    return X.t() * X / p;
}

arma::mat compute_k_matrix(const arma::mat& GSM,
                           const arma::rowvec& x_k,
                           int p) {
    return (GSM * p - arma::trans(x_k) * x_k) / (p - 1);
}

arma::mat compute_g_matrix(const arma::mat& K,
                           const arma::rowvec& x_k) {
    arma::mat G = K;
    G.each_row() %= x_k;
    G.each_col() %= arma::trans(x_k);
    return G;
}
