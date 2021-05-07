#!/bin/bash

cd "${MVMAPIT_DIR}" || exit # exit if 'cd ...' fails
git fetch --prune && git reset --hard origin/"${MVMAPIT_BRANCH_DEV}"
Rscript -e "install.packages('.', repos = NULL, type = 'source')"
