// Copyright 2021 Lorin Crawford.
#include <RcppArmadillo.h>
#include "mapit/davies.h"
#include "mqs/mqs.h"

// #define WITH_LOGGER 1 // uncomment for logging during development
#ifdef WITH_LOGGER  // check value
#define SPDLOG_DISABLE_DEFAULT_LOGGER 1
#include <RcppSpdlog>
#endif

arma::mat davies_routine(const arma::mat& S,
                         const arma::mat& Sinv,
                         const arma::mat& q,
                         const std::vector<arma::mat>& matrices) {
    arma::mat eigenvalues(matrices[0].n_cols, q.n_cols);

    for (int i = 0; i < q.n_cols; i++) {
        eigenvalues.col(i) = davies_routine_vec(S, Sinv, q.col(i), matrices);
    }
    return eigenvalues;
}

arma::vec davies_routine_vec(const arma::mat& S,
                             const arma::mat& Sinv,
                             const arma::vec& q,
                             const std::vector<arma::mat>& matrices) {
    int num_variance_components = matrices.size();
    arma::vec evals;
    arma::vec eigval;
    arma::mat eigvec;

    if (num_variance_components == 3) {  // C is NULL
        arma::vec q_sub = q.subvec(1, 2);
        arma::mat S_sub = S.submat(1, 1, 2, 2);
        arma::vec delta_null = arma::inv(S_sub) * q_sub;

        arma::eig_sym(eigval,
                      eigvec,
                      delta_null(0) * matrices[1]
                      + delta_null(1) * matrices[2]);
    } else {  // C is not NULL
        arma::vec delta_null = Sinv * q;

        arma::eig_sym(eigval,
                      eigvec,
                      delta_null(0) * matrices[1]
                      + delta_null(1) * matrices[2]
                      + delta_null(2) * matrices[3]);
    }
    arma::mat EV_mat = compute_positive_ev_matrix(eigvec, eigval);
    evals = arma::eig_sym(
                    EV_mat
                    * compute_h_matrix(Sinv, matrices)
                    * EV_mat);
    return evals;
}

arma::mat compute_positive_ev_matrix(const arma::mat& eigvec,
                                     const arma::mat& eigval) {
    arma::uvec ind_gt_zero = arma::find(eigval > 0);
    return eigvec.cols(ind_gt_zero)
            * arma::diagmat(sqrt(eigval(ind_gt_zero)))
            * arma::trans(eigvec.cols(ind_gt_zero));
}
