# mvMAPIT (development version)

* Added Cauchy combination test `cauchy_combined` including vignette that compares combination methods
* `simulate_traits` now returns genotype matrix with causal epistatic variants named according to the trait they affect
* Added progress bar and possibility to interrupt C++ routine using `RcppProgress`

# mvMAPIT 2.0.1 release

* Fix LTO issues when submitting to CRAN. The testthat issue https://github.com/r-lib/testthat/issues/1230
describes the solution chosen. Created this GitHub gist to reproduce LTO errors: https://gist.github.com/jdstamp/056475683110aacdb1e6761872ab1e05.

# mvMAPIT 2.0.0.1 pre-release

* CRAN issue does not show up in `R CMD check`. Version upgrade for resubmission.

# mvMAPIT 2.0.0 release

* `fishers_combined` and `harmonic_combined` now take additional arguments in 
form of string values. The first determines the name of the column of the tibble
by which to group the p-values. The second determines the name of the column 
  containing the p-values.
* Added `Dockerfile` and vignette on mvMAPIT in Docker.
* `MvMAPIT` was renamed to `mvmapit`.
* `simulate_phenotypes` was renamed to `simulate_traits`.
* Dependencies were cleaned up, e.g. R requirement 2.10 -> 3.5, old package 
dependencies removed.

# mvMAPIT 1.1.1-alpha prerelease
* Version that publication data was produced with. No official release.

