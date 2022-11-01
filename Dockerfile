FROM rocker/verse:4.0.5

WORKDIR ${HOME}
COPY ./ /tmp/mvMAPIT

RUN R -e "install.packages(c( 'checkmate', \
                              'CompQuadForm', \
                              'dplyr', \
                              'foreach', \
                              'harmonicmeanp', \
                              'logging', \
                              'mvtnorm', \
                              'Rcpp', \
                              'RcppAlgos', \
                              'RcppArmadillo', \
                              'RcppParallel', \
                              'RcppSpdlog', \
                              'stats', \
                              'testthat', \
                              'tidyr', \
                              'utils'), \
                             dependencies=TRUE \
                             );"

RUN R -e "install.packages('/tmp/mvMAPIT/', type = 'source', repos = NULL)"