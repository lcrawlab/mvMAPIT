# mvMAPIT (development version)


# mvMAPIT 2.0.0.1 release

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

