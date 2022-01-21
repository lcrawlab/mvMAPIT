// Copyright 2021 Lorin Crawford.
#include "mapit/pve.h"

arma::rowvec compute_pve(const arma::mat &variance_components,
                         int component_index) {
  arma::mat rowsum = arma::sum(variance_components, 0);
  return variance_components.row(component_index) / rowsum;
}
