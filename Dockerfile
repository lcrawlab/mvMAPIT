FROM rocker/verse:4.0.5

WORKDIR ${HOME}
COPY ./ /tmp/mvMAPIT
RUN mv /tmp/mvMAPIT/src/Makevars.in.dev /tmp/mvMAPIT/src/Makevars.in
RUN mv /tmp/mvMAPIT/src/MAPIT.h.dev /tmp/mvMAPIT/src/MAPIT.h

RUN git clone https://github.com/llvm-mirror/openmp.git
RUN cd openmp && mkdir build && cd build && cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ .. && make && make install
RUN apt-get install gcc g++ \
                      gfortran \
                      libblas-dev \
                      liblapack-dev \
                      libopenmpi-dev \
                      openmpi-bin

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

RUN R CMD INSTALL /tmp/mvMAPIT/ --preclean