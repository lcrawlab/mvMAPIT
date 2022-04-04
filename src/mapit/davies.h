// Copyright 2021 Lorin Crawford.
#pragma once
#include <RcppArmadillo.h>
#include <vector>

arma::mat davies_routine(const arma::mat &S, const arma::mat &Sinv,
                         const arma::mat &q,
                         const std::vector<arma::mat *> &matrices);

arma::vec davies_routine_vec(const arma::mat &S, const arma::mat &Sinv,
                             const arma::vec &q,
                             const std::vector<arma::mat *> &matrices);

arma::mat compute_positive_ev_matrix(const arma::mat &eigvec,
                                     const arma::mat &eigval);
