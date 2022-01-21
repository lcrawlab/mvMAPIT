// Copyright 2021 Lorin Crawford.
#include "logging/log.h"

std::string phenotype_vector_to_string(const std::vector<arma::vec> &v) {
  std::string s = vector_to_string(v[0]);
  for (int i = 1; i < v.size(); i++) {
    s.append("\n");
    s.append(vector_to_string(v[i]));
  }
  return s;
}

std::string vector_to_string(const arma::vec &v) {
  std::string s = std::to_string(v(0));
  for (int i = 1; i < v.n_elem; i++) {
    s.append(", ");
    s.append(std::to_string(v(i)));
  }
  return s;
}

std::string matrix_to_string(const arma::mat &matrix) {
  int n = matrix.n_rows;
  int m = matrix.n_cols;
  std::string s = "";
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      s.append(std::to_string(matrix(i, j)));
      s.append(", ");
    }
    s.append("\n");
  }
  return s;
}
