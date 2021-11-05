// Copyright 2021 Lorin Crawford.
#include "mapit/normal.h"

// #define WITH_LOGGER 1  // uncomment for logging during development

#ifdef WITH_LOGGER  // check value
#define SPDLOG_DISABLE_DEFAULT_LOGGER 1
#include <RcppSpdlog>
#include <spdlog/stopwatch.h>  // also support stopwatch feature
#endif

arma::mat normal_pvalues(const arma::mat& variance_estimate,
                                   const arma::mat& standard_error) {
#ifdef WITH_LOGGER
    std::string logname = "mapit.normal.normal_pvalues";
    auto logger = spdlog::get(logname);  // retrieve existing one
     // or create new one if needed
    if (logger == nullptr) logger = spdlog::r_sink_mt(logname);
#endif
    arma::mat pvalues(variance_estimate.n_rows, variance_estimate.n_cols);
    arma::mat ratio(variance_estimate.n_rows, variance_estimate.n_cols);
    ratio = abs(variance_estimate / standard_error);
    for (int j = 0; j < variance_estimate.n_cols; j++) {
        for (int i = 0; i < variance_estimate.n_rows; i++) {
            pvalues(i, j) =
                2 * R::pnorm(arma::as_scalar(ratio(i, j)), 0.0, 1.0, 0, 0);
#ifdef WITH_LOGGER
                logger->info("ratio({}, {}) = {}", i, j, ratio(i, j));
                logger->info("pvalues({}, {}) = {}", i, j, pvalues(i, j));
                logger->info("var({}, {}) = {}", i, j, variance_estimate(i, j));
                logger->info("se({}, {}) = {}", i, j, standard_error(i, j));
#endif
        }
    }
    return pvalues;
}
