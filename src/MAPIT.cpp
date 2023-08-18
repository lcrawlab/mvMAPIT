// Copyright 2017-2021 Lorin Crawford.

#include "MAPIT.h"
#include "gsm/gsm.h"
#include "logging/log.h"
#include "mapit/davies.h"
#include "mapit/normal.h"
#include "mapit/projection.h"
#include "mapit/pve.h"
#include "mapit/util.h"
#include "mqs/mqs.h"

using std::chrono::duration_cast;

using std::chrono::milliseconds;
using std::chrono::steady_clock;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

////////////////////////////////////////////////////////////////////////////

// Below are functions for MAPIT using two hypothesis testing strategies:
// (1) MAPIT using the Normal or Z-Test
// (2) MAPIT using the Davies Method

// Considered are the following submodels:
// (1) Standard Model ---> y = m+g+e
// (2) Standard + Covariate Model ---> y = Wa+m+g+e
// (3) Standard + Common Environment Model ---> y = m+g+c+e
// (4) Standard + Covariate + Common Environment Model ---> y = Wa+m+g+c+e

// NOTE: delta = {delta(0),delta(1),delta(2)} = {sigma^2,omega^2,tau^2}
// for models (1) and (2)
// NOTE: delta = {delta(0),delta(1),delta(2),delta(3)}
//             = {sigma^2,omega^2,nu^2,tau^2} for models (3) and (4)

////////////////////////////////////////////////////////////////////////////

// Generalized MAPIT -- should handle all models

// [[Rcpp::export]]
Rcpp::List MAPITCpp(
        const arma::mat X, const arma::mat Y,
        Rcpp::Nullable <Rcpp::NumericMatrix> Z = R_NilValue,
        Rcpp::Nullable <Rcpp::NumericMatrix> C = R_NilValue,
        Rcpp::Nullable <Rcpp::NumericVector> variantIndices = R_NilValue,
        std::string testMethod = "normal", int cores = 1,
        Rcpp::Nullable <Rcpp::NumericMatrix> GeneticSimilarityMatrix = R_NilValue) {
    int i;
    const int n = X.n_cols;
    const int p = X.n_rows;
    const int d = Y.n_rows;
    int num_combinations = 1;
    int z = 0;

    const bool combinatorial = (d > 1);
    if (combinatorial) {
        num_combinations = num_combinations_with_replacement(d, 2);
    }

#ifdef WITH_LOGGER
    std::string logname = "MAPITcpp";
    auto logger = spdlog::get(logname);
    if (logger == nullptr)
      logger = spdlog::r_sink_mt(logname);

    // instantiate a stop watch
    spdlog::stopwatch sw;

    logger->info("Number of samples: {}", n);
    logger->info("Number of SNPs: {}", p);
    logger->info("Number of phenotypes: {}", d);
    logger->info("Test method: {}", testMethod);

#ifdef _OPENMP
    logger->info("Execute c++ routine on {} cores.", cores);
#endif

#endif

    arma::mat sigma_est(p, num_combinations);
    arma::mat sigma_se(p, num_combinations);
    arma::mat pve(p, num_combinations);
    arma::mat execution_t(p, 6);

    int L_rows, L_cols;
    if (combinatorial) {
        L_rows = n;
        L_cols = num_combinations;
    } else {
        L_rows = n * d;
        L_cols = 1;
    }
    arma::cube Lambda(L_rows, L_cols, p);

    arma::mat GSM;
    if (GeneticSimilarityMatrix.isNull()) {
        GSM = get_linear_kernel(X);
    } else {
        GSM = Rcpp::as<arma::mat>(GeneticSimilarityMatrix.get());
    }

    arma::mat Zz;
    if (Z.isNotNull()) {
        Zz = Rcpp::as<arma::mat>(Z.get());
        z = Zz.n_rows;
    }

    arma::vec ind;
    // check if we are provided variants of interest
    if (variantIndices.isNotNull()) {
        ind = Rcpp::as<arma::vec>(variantIndices.get());
    }

#ifdef _OPENMP
    omp_set_num_threads(cores);
#endif
    Progress prgrs(p, true);
#pragma omp parallel for schedule(dynamic)
    for (i = 0; i < p; i++) {
        if (!Progress::check_abort()) {
            prgrs.increment();
#ifdef WITH_LOGGER
            logger->info("Variant {}/{}.", i + 1, p);
#endif

            if (skip_variant(ind, i)) {
#ifdef WITH_LOGGER
                logger->info("Variant {}/{} not of interest. Skip to next.", i + 1, p);
#endif
                continue;
            }

            std::vector <arma::mat> matrices;
            matrices.resize(4);
            arma::mat &G = matrices[0];
            arma::mat &K = matrices[1];
            arma::mat &M = matrices[2];
            arma::mat &Cc = matrices[3];
            if (C.isNotNull()) {
                Cc = Rcpp::as<arma::mat>(C.get());
            } else {
                matrices.resize(3);
            }

            // Compute K and G covariance matrices
            auto start = steady_clock::now();
            // Create the linear kernel
            const arma::rowvec x_k = X(arma::span(i), arma::span::all);
            K = compute_k_matrix(GSM, x_k, p);
            G = compute_g_matrix(K, x_k);

            auto end = steady_clock::now();
            execution_t(i, 0) = duration_cast<milliseconds>(end - start).count();

#ifdef WITH_LOGGER_FINE
            logger->info("Dimensions of polygenic background: {} x {}.", K.n_cols,
                         K.n_rows);
#endif

            // Transform K and G using projection M
            start = steady_clock::now();
            arma::mat b = arma::zeros(n, z + 2);
            b.col(0) = arma::ones<arma::vec>(n);
            if (z > 0) {
                b.cols(1, z) = Zz.t();
            }
            b.col(z + 1) = arma::trans(x_k);

            M = compute_projection_matrix(n, b);
            K = project_matrix(K, b);
            G = project_matrix(G, b);

            if (C.isNotNull()) {
                Cc = project_matrix(Cc, b);
            }
            b.reset();

            const arma::mat Yc = Y * M;
            end = steady_clock::now();
            execution_t(i, 1) = duration_cast<milliseconds>(end - start).count();

            arma::mat q;
            std::vector <arma::vec> phenotypes;
            start = steady_clock::now();
            if (combinatorial) {
                phenotypes = matrix_to_vector_of_rows(Yc);

            } else {
                arma::vec yc = vectorise(Yc);
                phenotypes = matrix_to_vector_of_rows(yc.as_row());
            }
            end = steady_clock::now();
            execution_t(i, 2) = duration_cast<milliseconds>(end - start).count();

            start = steady_clock::now();
            q = compute_q_matrix(phenotypes, matrices);
            end = steady_clock::now();
            execution_t(i, 3) = duration_cast<milliseconds>(end - start).count();

            start = steady_clock::now();
            arma::mat S = compute_s_matrix(matrices);
            end = steady_clock::now();
            execution_t(i, 4) = duration_cast<milliseconds>(end - start).count();

            // Compute delta and Sinv
            const float det_S = arma::det(S);
            if (det_S == 0) {
#ifdef WITH_LOGGER
                logger->warn("The determinant of the S matrix is {}.", det_S);
                logger->info("Skip variant {}.", i + 1);
#endif
                continue;
            }
            arma::mat Sinv = arma::inv(S);
            arma::mat delta = Sinv * q;

            // Save point estimates of the epistasis component
            sigma_est.row(i) = delta.row(0);
#ifdef WITH_LOGGER_FINE
            logger->info("S-1({}):\n {}.", i + 1, matrix_to_string(Sinv));
            logger->info("delta({}):\n {}.", i + 1, matrix_to_string(delta));
            logger->info("q({}):\n {}.", i + 1, matrix_to_string(q));
#endif

            start = steady_clock::now();
            if (testMethod == "normal") {
                arma::vec var_delta =
                        compute_variance_delta(phenotypes, Sinv, delta, matrices);
                // Save SE of the epistasis component
                sigma_se.row(i) = arma::trans(sqrt(var_delta));
#ifdef WITH_LOGGER_FINE
                logger->info("var_delta({}):\n {}.", i + 1,
                             vector_to_string(sqrt(var_delta)));
#endif
            } else if (testMethod == "davies") {
                try {
                    Lambda.slice(i) = davies_routine(S, Sinv, q, matrices);

#ifdef WITH_LOGGER_FINE
                    logger->info("Lambda.slice({}). {}", i + 1,
                                 matrix_to_string(Lambda.slice(i)));
#endif
                } catch (std::exception &e) {
#ifdef WITH_LOGGER
                    logger->error("Error: {}.", e.what());
                    logger->info("Skip davies method for variant {}.", i + 1);
#endif
                    continue;
                }
            }
            end = ::steady_clock::now();
            execution_t(i, 5) = duration_cast<milliseconds>(end - start).count();
            // Compute the PVE
            pve.row(i) = compute_pve(delta, 0);
#ifdef WITH_LOGGER_FINE
            logger->info("PVE({}): \n{}.", i + 1, matrix_to_string(pve.row(i)));
#endif
        }
    }
#ifdef WITH_LOGGER
    logger->info("Elapsed time: {}", sw);
#endif
    // Return a Rcpp::List of the arguments
    if (testMethod == "davies") {
#ifdef WITH_LOGGER
        logger->info("Return from davies method.");
#endif
        return Rcpp::List::create(
                Rcpp::Named("Est") = sigma_est, Rcpp::Named("Eigenvalues") = Lambda,
                Rcpp::Named("PVE") = pve, Rcpp::Named("timings") = execution_t);
    } else {
#ifdef WITH_LOGGER
        logger->info("Return from normal method.");
#endif
        // Compute the p-values for each estimate
        arma::mat pvalues = normal_pvalues(sigma_est, sigma_se);
        // H0: sigma = 0 vs. H1: sigma != 0
#ifdef WITH_LOGGER_FINE
        logger->info("sigma_est({}):\n {}.", i + 1, matrix_to_string(sigma_est));
        logger->info("sigma_se({}):\n {}.", i + 1, matrix_to_string(sigma_se));
        logger->info("pvalues({}):\n {}.", i + 1, matrix_to_string(pvalues));
#endif
        return Rcpp::List::create(
                Rcpp::Named("Est") = sigma_est, Rcpp::Named("SE") = sigma_se,
                Rcpp::Named("pvalues") = pvalues, Rcpp::Named("PVE") = pve,
                Rcpp::Named("timings") = execution_t);
    }
}
