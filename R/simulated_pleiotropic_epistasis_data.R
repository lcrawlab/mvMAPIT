#' Genotype and multiple phenotype data with pleiotropic epistasis.
#'
#' A simulated dataset that has pleiotropic epistatic interactions.
#'
#' @format A list with simulated data and simulation parameters:
#' \describe{
#'   \item{numer_snp}{Number of SNPs for each genotype.}
#'   \item{number_samples}{Number of genotype samples.}
#'   \item{pve}{Phenotypic variance explained/broad-sense heritability (H^2).}
#'   \item{rho}{Portion of H^2 that is contributed to by the marginal (additive) effects.}
#'   \item{phenotype}{The simulated phenotype. Vector size (number_samples) x 1.}
#'   \item{genotype}{File name of the genotype data.}
#'   \item{causal_snps}{Nested list of causal SNP IDs and effect sizes.}
#' }
#' @source data-raw/simulate_epistasis.R
"simulated_pleiotropic_epistasis_data"
