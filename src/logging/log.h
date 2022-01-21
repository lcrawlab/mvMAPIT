// Copyright 2021 Lorin Crawford.
#pragma once
#include <RcppArmadillo.h>
#include <string>
#include <vector>

std::string phenotype_vector_to_string(const std::vector<arma::vec> &v);

std::string vector_to_string(const arma::vec &v);

std::string matrix_to_string(const arma::mat &matrix);
