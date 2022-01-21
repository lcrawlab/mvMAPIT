// Copyright 2021 Lorin Crawford.
#pragma once
#include <RcppArmadillo.h>

arma::rowvec compute_pve(const arma::mat &variance_components,
                         int component_index);
