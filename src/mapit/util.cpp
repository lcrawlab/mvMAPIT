// Copyright 2021 Lorin Crawford.
#include <RcppArmadillo.h>
#include "mapit/util.h"

arma::mat compute_principal_components(const arma::mat& X, int top = 10) {
    arma::mat U;
    arma::vec s;
    arma::mat V;
    svd(U, s, V, X);

    arma::mat PCs = U * diagmat(s);
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

std::vector<arma::vec> matrix_to_vector_of_rows(const arma::mat& matrix) {
    std::vector<arma::vec> V(matrix.n_rows);
    for (int i = 0; i < matrix.n_rows; ++i) {
        V[i] = arma::conv_to<arma::vec>::from(matrix.row(i));
    }
    return V;
}

arma::mat vectorise_to_matrix(const arma::mat& matrix) {
    arma::mat V(matrix.n_rows * matrix.n_cols, 1);
    V.col(0) = vectorise(matrix);
    return V;
}

int factorial(int n) {
     return (n == 0) || (n == 1) ? 1 : n * factorial(n-1);
}

int num_combinations_with_replacement(int num_available, int num_selected) {
     return factorial(num_available + num_selected - 1)
                / (factorial(num_selected) * factorial(num_available - 1));
}
