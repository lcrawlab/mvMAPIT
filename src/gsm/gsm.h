// Copyright 2021 Lorin Crawford.
#pragma once

arma::mat GetLinearKernel(const arma::mat X);

arma::mat ComputeKMatrix(const arma::mat& GSM, const arma::rowvec& x_k, int p);

arma::mat ComputeGMatrix(const arma::mat& K, const arma::rowvec& x_k);
