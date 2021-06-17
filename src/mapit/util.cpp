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
