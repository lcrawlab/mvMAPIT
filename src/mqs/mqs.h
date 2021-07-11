// Copyright 2021 Lorin Crawford.
#pragma once
#include <RcppArmadillo.h>
#include <vector>
#include <string>

double product_trace(const arma::mat& a, const arma::mat& b);

arma::mat compute_s_matrix(const std::vector<arma::mat>& matrices);

arma::vec compute_q_vector(const arma::vec& y,
                           const std::vector<arma::mat>& matrices);

arma::vec compute_q_vector(const arma::vec& y1,
                           const arma::vec& y2,
                           const std::vector<arma::mat>& matrices);

arma::mat compute_q_matrix(const std::vector<arma::vec>& Y,
                           const std::vector<arma::mat>& matrices);

arma::mat compute_h_matrix(const arma::mat& Sinv,
                           const std::vector<arma::mat>& matrices);

arma::mat compute_v_matrix(const arma::vec& delta,
                           const std::vector<arma::mat>& matrices);

double compute_mqs_var_approximation(const arma::vec& yc,
                              const arma::mat& H,
                              const arma::mat& V);

double compute_mqs_var_approximation(const arma::vec& y1,
                              const arma::vec& y2,
                              const arma::mat& H,
                              const arma::mat& V);

double compute_variance_delta(const arma::vec& yc,
                              const arma::mat& Sinv,
                              const arma::vec& delta,
                              const std::vector<arma::mat>& matrices);

arma::vec compute_variance_delta(const std::vector<arma::vec>& Y,
                                 const arma::mat& Sinv,
                                 const arma::mat& delta,
                                 const std::vector<arma::mat>& matrices);
