// Copyright 2021 Lorin Crawford.
#include "mapit/kronecker.h"

std::vector<arma::mat> kronecker_products(std::vector<arma::mat> matrices,
                                         const arma::mat V_phenotype,
                                         const arma::mat V_error) {
    int num_variance_components = matrices.size();
    for (int i = 0; i < num_variance_components - 1; i++) {
        matrices[i] = kron(V_phenotype, matrices[i]);
    }
    matrices[num_variance_components - 1]
            = kron(V_error, matrices[num_variance_components - 1]);
    return matrices;
}
