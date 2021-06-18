// Copyright 2021 Lorin Crawford.
#pragma once

arma::vec davies_routine(const arma::mat& S,
                             const arma::mat& Sinv,
                             const arma::vec& q,
                             const std::vector<arma::mat>& matrices);

arma::mat compute_positive_ev_matrix(const arma::mat& eigvec,
                                     const arma::mat& eigval);
