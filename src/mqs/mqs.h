// Copyright 2021 Lorin Crawford.
#pragma once
#include <vector>
#include <string>

double product_trace(const arma::mat& a, const arma::mat& b);

arma::mat compute_s_matrix(const std::vector<arma::mat>& matrices);

arma::vec compute_q_vector(const arma::vec& yc,
                           const std::vector<arma::mat>& matrices);

arma::mat compute_h_matrix(const arma::mat& Sinv,
                           const std::vector<arma::mat>& matrices);

arma::mat compute_v_matrix(const arma::vec& delta,
                           const std::vector<arma::mat>& matrices);

double compute_variance_delta(const arma::vec& yc,
                              const arma::mat& H,
                              const arma::mat& V);

double compute_variance_delta(const arma::vec& yc,
                              const arma::mat& Sinv,
                              const arma::vec& delta,
                              const std::vector<arma::mat>& matrices);
