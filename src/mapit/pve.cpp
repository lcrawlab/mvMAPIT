// Copyright 2021 Lorin Crawford.
#include "mapit/pve.h"

arma::vec compute_pve(arma::mat variance_components) {
    arma::rowvec rowsum = arma::sum(variance_components, 0);
    return variance_components.row(0) / rowsum;
}

