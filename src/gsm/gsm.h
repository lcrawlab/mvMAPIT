// Copyright 2021 Lorin Crawford.
#pragma once

arma::mat get_linear_kernel(const arma::mat& X);

arma::mat compute_k_matrix(const arma::mat& GSM,
                           const arma::rowvec& x_k,
                           int p);

arma::mat compute_g_matrix(const arma::mat& K,
                           const arma::rowvec& x_k);
