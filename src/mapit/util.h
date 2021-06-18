// Copyright 2021 Lorin Crawford.
#pragma once

arma::mat compute_principal_components(const arma::mat& X, int top);

arma::vec remove_first_element(const arma::vec& vector);

bool skip_variant(const arma::vec& ind, int i);

std::vector< arma::vec > matrix_to_vector_of_vectors(const arma::mat& matrix);
