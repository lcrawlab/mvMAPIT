// Copyright 2021 Lorin Crawford.
#include <RcppArmadillo.h>
#include "mapit/util.h"

arma::mat compute_principal_components(const arma::mat& X, int top = 10) {
    arma::mat U;
    arma::vec s;
    arma::mat V;
    svd(U, s, V, X);

    arma::mat PCs = U*diagmat(s);
    return PCs.cols(0, top - 1);
}

arma::vec remove_first_element(const arma::vec& vector) {
    arma::vec new_vector = arma::vec(vector.size() - 1);
    for (int i = 1; i < vector.size(); i++) {
        new_vector(i - 1) = vector(i);
    }
    return new_vector;
}

bool skip_variant(const arma::vec& ind, int i) {
    if (ind.is_empty()) {
        // look for i+1 because R uses 1-based indexing
        // if there is no match find == ind.end()
        return false;
    }
    if (std::find(ind.begin(), ind.end(), i+1) == ind.end()) {
        return true;
    }
    return false;
}
