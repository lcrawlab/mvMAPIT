#ifndef MAPIT_H
#define MAPIT_H

#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#define ARMA_64BIT_WORD 1
#include <vector>
#include <string>
#include <iostream>

#include <Rcpp.h>
using namespace Rcpp;

arma::mat GetLinearKernel(arma::mat X);
arma::mat ComputePCs(arma::mat X,int top);
//TODO: tests
arma::mat ComputeProjectionMatrix(int n, arma::mat b);
double ProductTrace(arma::mat a, arma::mat b);

#endif