FROM rocker/verse:4.0.5

COPY ./ /tmp/mvMAPIT

RUN R -e "install.packages('/tmp/mvMAPIT/', type = 'source', repos = NULL, dependencies = TRUE)"
