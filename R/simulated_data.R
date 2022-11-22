#' Genotype and trait data with epistasis.
#'
#' A simulated dataset that has epistatic interactions.
#'
#' @format A named list with simulated data and simulation parameters:
#' \describe{
#'   \item{parameters}{Tibble containing simulation parameters for each trait.}
#'   \item{trait}{Matrix containing simulated data for 2 traits and 500 samples.}
#'   \item{genotype}{Matrix containing simulated genotype with 500 samples and 1000 variables.}
#'   \item{additive}{Tibble containing all variants with additive effects on the traits as well as the effect sizes.}
#'   \item{epistatic}{Tibble containing all variants with epistatic effects on the traits as well as the effect sizes.}
#'   \item{interactions}{Tibble containing all interactions, effect size, and trait they affect.}
#'   \item{snps.filtered}{SNPs that were used in the simulations.}
#' }
#' @source data-raw/simulate_epistasis.R
"simulated_data"
