// Copyright 2021 Lorin Crawford.
#pragma once

std::vector<arma::mat> kronecker_products(std::vector<arma::mat> matrices,
                                         const arma::mat V_phenotype,
                                         const arma::mat V_error);
