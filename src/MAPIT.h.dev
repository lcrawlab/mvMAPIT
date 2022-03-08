// Copyright 2017-2021 Lorin Crawford.
#pragma once

#include <exception>
#include <string>
#include <typeinfo>
#include <vector>

#define WITH_LOGGER 1  // uncomment for logging during development
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef WITH_LOGGER
#define SPDLOG_DISABLE_DEFAULT_LOGGER 1
#include <RcppSpdlog>
#include <spdlog/stopwatch.h>
#endif
