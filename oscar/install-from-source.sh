#!/bin/bash

module load R/4.0.5
module load gcc/10.2 pcre2/10.35 intel/2020.2 texlive/2018
module load lapack/3.7.0 openblas/0.3.7

cd "${MVMAPIT_DIR}" || exit # exit if 'cd ...' fails
git fetch --prune && git reset --hard origin/"${MVMAPIT_BRANCH_DEV}"
Rscript -e "install.packages('.', repos = NULL, type = 'source')"
