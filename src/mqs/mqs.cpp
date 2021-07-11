// Copyright 2021 Lorin Crawford.
#include "mqs/mqs.h"
// #define WITH_LOGGER 1  // uncomment for logging during development
#ifdef WITH_LOGGER  // check value
#define SPDLOG_DISABLE_DEFAULT_LOGGER 1
#include <RcppSpdlog>
#endif
#include "mapit/util.h"

// Computes Tr(ab)
double product_trace(const arma::mat& a, const arma::mat& b) {
    // Computed efficiently using Hadamard (elementwise) product
    // https://en.wikipedia.org/wiki/Trace_(linear_algebra)#Trace_of_a_product
    // https://proofwiki.org/wiki/Trace_of_Matrix_Product
    return arma::as_scalar(accu(a%b));
}

arma::mat compute_s_matrix(const std::vector<arma::mat>& matrices) {
#ifdef WITH_LOGGER
    std::string logname = "mqs::mqs::compute_s_matrix";
    auto logger = spdlog::get(logname);  // retrieve existing one
     // or create new one if needed
    if (logger == nullptr) logger = spdlog::r_sink_mt(logname);
#endif
    int num_variance_components = matrices.size();
    arma::mat S = arma::zeros(num_variance_components, num_variance_components);

    for (int i = 0; i < num_variance_components; i++) {
        for (int j = 0; j < num_variance_components; j++) {
            if (i <= j) {  // create upper triangular matrix
                S(i, j) = product_trace(matrices[i], matrices[j]);
                S(j, i) = S(i, j);
#ifdef WITH_LOGGER
                logger->info("S({},{}) = {}", i, j, S(i, j));
                if (S(i, j) < 0) {
                    logger->warn("S({},{}) negative", i, j);
                }
#endif
            }
        }
    }
    return S;
}

arma::vec compute_q_vector(const arma::vec& y,
                           const std::vector<arma::mat>& matrices) {
    return compute_q_vector(y, y, matrices);
}

arma::vec compute_q_vector(const arma::vec& y1,
                           const arma::vec& y2,
                           const std::vector<arma::mat>& matrices) {
#ifdef WITH_LOGGER
    std::string logname = "mqs::mqs::compute_q_vector";
    auto logger = spdlog::get(logname);
    if (logger == nullptr) logger = spdlog::r_sink_mt(logname);
#endif
    int num_variance_components = matrices.size();
    arma::vec q = arma::zeros(num_variance_components);

    for (int i = 0; i < num_variance_components; i++) {
        q(i) = arma::as_scalar(y1.t() * matrices[i] * y2);
#ifdef WITH_LOGGER
        logger->info("q({}) = {}", i, q(i));
        if (q(i) < 0) {
            logger->warn("q({}) negative", i);
        }
#endif
    }
    return q;
}

arma::mat compute_q_matrix(const std::vector<arma::vec>& Y,
                           const std::vector<arma::mat>& matrices) {
#ifdef WITH_LOGGER
    std::string logname = "mqs::mqs::compute_q_matrix";
    auto logger = spdlog::get(logname);
    if (logger == nullptr) logger = spdlog::r_sink_mt(logname);
    logger->info("Computing q matrix");
#endif
    int num_variance_components = matrices.size();
    int max_index = Y.size();
    int num_combinations =
        num_combinations_with_replacement(max_index, 2);
    int index_q_col = 0;
    arma::mat q; q.zeros(num_variance_components, num_combinations);
#ifdef WITH_LOGGER
    logger->info("Number variance components: {}", num_variance_components);
    logger->info("Max index Y vector: {}", max_index);
    logger->info("Number of combinations: {}", num_combinations);
#endif
    for (int j = 0; j < max_index; j++) {
        for (int k = 0; k < max_index; k++) {
            if (k <= j) {
    #ifdef WITH_LOGGER
                logger->info("Combination: {}, {}", j, k);
                logger->info("q column: {}", index_q_col);
    #endif
                q.col(index_q_col) = compute_q_vector(Y[j], Y[k], matrices);
                index_q_col += 1;
            }
        }
    }
    return q;
}

arma::mat compute_h_matrix(const arma::mat& Sinv,
                           const std::vector<arma::mat>& matrices) {
#ifdef WITH_LOGGER
    std::string logname = "mqs::mqs::compute_h_matrix";
    auto logger = spdlog::get(logname);
    if (logger == nullptr) logger = spdlog::r_sink_mt(logname);
    logger->info("Computing H matrix");
#endif
    arma::mat H = arma::zeros(arma::size(matrices[0]));
    for (int i = 0; i < matrices.size(); i++) {
        H = H + Sinv(0, i) * matrices[i];
    }
    return H;
}

arma::mat compute_v_matrix(const arma::vec& delta,
                           const std::vector<arma::mat>& matrices) {
#ifdef WITH_LOGGER
    std::string logname = "mqs::mqs::compute_v_matrix";
    auto logger = spdlog::get(logname);
    if (logger == nullptr) logger = spdlog::r_sink_mt(logname);
    logger->info("Computing V matrix");
#endif
    arma::mat V = arma::zeros(arma::size(matrices[0]));
    for (int i = 0; i < matrices.size(); i++) {
        V = V + delta(i) * matrices[i];
    }
    return V;
}

double compute_mqs_var_approximation(const arma::vec& yc,
                              const arma::mat& H,
                              const arma::mat& V) {
#ifdef WITH_LOGGER
    std::string logname = "mqs::mqs::compute_mqs_var_approximation";
    auto logger = spdlog::get(logname);
    if (logger == nullptr) logger = spdlog::r_sink_mt(logname);
    logger->info("Computing variance of variance component. One phenotype.");
#endif
    arma::mat Hy = H * yc;
    return arma::as_scalar(2 * Hy.t() * V * Hy);
}

double compute_mqs_var_approximation(const arma::vec& y1,
                              const arma::vec& y2,
                              const arma::mat& H,
                              const arma::mat& V) {
    return arma::as_scalar(2 * y1.t() * H.t() * V * H * y2);
}

double compute_variance_delta(const arma::vec& yc,
                              const arma::mat& Sinv,
                              const arma::vec& delta,
                              const std::vector<arma::mat>& matrices) {
    arma::mat H = compute_h_matrix(Sinv, matrices);
    arma::mat V = compute_v_matrix(delta, matrices);
#ifdef WITH_LOGGER
    std::string logname = "mqs::mqs::compute_variance_delta";
    auto logger = spdlog::get(logname);
    if (logger == nullptr) logger = spdlog::r_sink_mt(logname);
    const float det_H = arma::det(H);
    if (det_H == 0) {
        logger->warn("The determinant of the H matrix is {}.", det_H);
    }
    const float det_V = arma::det(V);
    if (det_V == 0) {
        logger->warn("The determinant of the V matrix is {}.", det_V);
    }
#endif
    return compute_mqs_var_approximation(yc, H, V);
}

arma::vec compute_variance_delta(const std::vector<arma::vec>& Y,
                                 const arma::mat& Sinv,
                                 const arma::mat& delta,
                                 const std::vector<arma::mat>& matrices) {
#ifdef WITH_LOGGER
    std::string logname = "mqs.mqs.compute_variance_delta";
    auto logger = spdlog::get(logname);
    if (logger == nullptr) logger = spdlog::r_sink_mt(logname);
    logger->info("Computing variance of variance component. Two phenotypes.");
#endif
    arma::mat H = compute_h_matrix(Sinv, matrices);
    arma::vec variance(delta.n_cols, arma::fill::zeros);
    int max_index = Y.size();
    int index_delta = 0;
    for (int j = 0; j < max_index; j++) {
        for (int k = 0; k < max_index; k++) {
            if (k <= j) {
#ifdef WITH_LOGGER
            logger->info("({},{}) {}/{}.", j,
                                           k,
                                           index_delta + 1,
                                           delta.n_cols);
#endif
                arma::mat V = compute_v_matrix(delta.col(index_delta),
                                                matrices);
                variance(index_delta) = compute_mqs_var_approximation(Y[j],
                                                                      Y[k],
                                                                      H,
                                                                      V);
                index_delta += 1;
            }
        }
    }
    return variance;
}
