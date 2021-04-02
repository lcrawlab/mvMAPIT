#' Return values of MAPIT applied to simulated_epistasis_data.rda.
#'
#' A dataset containing the p-values of different MAPIT runs.
#'
#' @format A list with p-values and interaction SNPs:
#' \describe{
#'   \item{davies_pvalues}{The p-values returned by MAPIT using the "davies" method.}
#'   \item{hybrid_pvalues}{The p-values returned by MAPIT using the "hybrid" method.}
#'   \item{normal_pvalues}{The p-values returned by MAPIT using the "normal" method.}
#'   \item{exhaustive_search}{A list containing the p-values of an exhaustive search and fitting the interaction of significant SNPs.}
#'   \item{epistatic_snps}{The list of the epistatic SNPs.}
#'   \item{additive_snps}{The list of the additive SNPs.}
#' }
#' @source data-raw/mapit_on_simulated_data.R
"mapit_analysis_data"
