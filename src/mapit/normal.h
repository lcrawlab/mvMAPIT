// Copyright 2021 Lorin Crawford.
#pragma once
#include <RcppArmadillo.h>
#include <string>

arma::mat normal_pvalues(const arma::mat &variance_estimate,
                         const arma::mat &standard_error);
