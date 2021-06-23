// Copyright 2021 Lorin Crawford.
#pragma once
#include <RcppArmadillo.h>
#include <testthat.h>
#include <typeinfo>
#include <vector>

arma::mat compute_principal_components(const arma::mat& X, int top);

arma::vec remove_first_element(const arma::vec& vector);

bool skip_variant(const arma::vec& ind, int i);

std::vector< arma::vec > matrix_to_vector_of_rows(const arma::mat& matrix);

arma::mat vectorise_to_matrix(const arma::mat& matrix);

int factorial(int n);

int num_combinations_with_replacement(int num_available, int num_selected);
