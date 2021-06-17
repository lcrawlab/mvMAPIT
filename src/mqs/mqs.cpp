// Copyright 2021 Lorin Crawford.
#include <RcppArmadillo.h>
#include "mqs.h"

// #define WITH_LOGGER 1 // uncomment for logging during development
#ifdef WITH_LOGGER  // check value
#define SPDLOG_DISABLE_DEFAULT_LOGGER 1
#include <RcppSpdlog>
#endif

// Computes Tr(ab)
double ProductTrace(arma::mat a, arma::mat b) {
    // Computed efficiently using Hadamard (elementwise) product
    // https://en.wikipedia.org/wiki/Trace_(linear_algebra)#Trace_of_a_product
    // https://proofwiki.org/wiki/Trace_of_Matrix_Product
    return arma::as_scalar(accu(a%b));
}

arma::mat ComputeSMatrix(std::vector<arma::mat> matrices) {
#ifdef WITH_LOGGER
    std::string logname = "mqs::mqs::ComputeSMatrix";
    auto logger = spdlog::get(logname);  // retrieve existing one
     // or create new one if needed
    if (logger == nullptr) logger = spdlog::r_sink_mt(logname);
#endif
    int num_variance_components = matrices.size();
    arma::mat S = arma::zeros(num_variance_components, num_variance_components);

    for (int i = 0; i < num_variance_components; i++) {
        for (int j = 0; j < num_variance_components; j++) {
            if (i <= j) {  // create upper triangular matrix
                S(i, j) = ProductTrace(matrices[i], matrices[j]);
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

arma::vec ComputeqVector(arma::vec yc, std::vector<arma::mat> matrices) {
#ifdef WITH_LOGGER
    std::string logname = "mqs::mqs::ComputeqVector";
    auto logger = spdlog::get(logname);
    if (logger == nullptr) logger = spdlog::r_sink_mt(logname);
#endif
    int num_variance_components = matrices.size();
    arma::vec q = arma::zeros(num_variance_components);

    for (int i = 0; i < num_variance_components; i++) {
        q(i) = arma::as_scalar(yc.t() * matrices[i] * yc);
#ifdef WITH_LOGGER
        logger->info("q({}) = {}", i, q(i));
        if (q(i) < 0) {
            logger->warn("q({}) negative", i);
        }
#endif
    }
    return q;
}

// TODO(jdstamp) ComputeHMatrix and ComputeVMatrix can probably be merged
arma::mat ComputeHMatrix(arma::mat Sinv, std::vector<arma::mat> matrices) {
#ifdef WITH_LOGGER
    std::string logname = "mqs::mqs::ComputeHMatrix";
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

arma::mat ComputeVMatrix(arma::vec delta, std::vector<arma::mat> matrices) {
#ifdef WITH_LOGGER
    std::string logname = "mqs::mqs::ComputeVMatrix";
    auto logger = spdlog::get(logname);
    if (logger == nullptr) logger = spdlog::r_sink_mt(logname);
    logger->info("Computing V matrix");
#endif
    arma::mat V = arma::zeros(arma::size(matrices[0]));
    for (int i = 0; i < matrices.size(); i++) {
        V = V + delta(i)*matrices[i];
    }
    return V;
}
