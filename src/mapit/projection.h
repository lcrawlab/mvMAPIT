// Copyright 2021 Lorin Crawford.
#pragma once
#include <RcppArmadillo.h>

arma::mat compute_projection_matrix(int n, const arma::mat &b);

arma::mat project_matrix(const arma::mat &A, const arma::mat &b);
