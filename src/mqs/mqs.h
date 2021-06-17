// Copyright 2021 Lorin Crawford.
#pragma once
double ProductTrace(arma::mat a, arma::mat b);

arma::mat ComputeSMatrix(std::vector<arma::mat> matrices);

arma::vec ComputeqVector(arma::vec yc, std::vector<arma::mat> matrices);

arma::mat ComputeHMatrix(arma::mat Sinv, std::vector<arma::mat> matrices);

arma::mat ComputeVMatrix(arma::vec delta, std::vector<arma::mat> matrices);
