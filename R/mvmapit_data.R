#' Multivariate MAPIT analysis and exhaustive search analysis.
#'
#' This data set contains the return object from the multivariate MAPIT method,
#' the fisher combined p-values, and the result from an exhaustive search using
#' regression on the SNPs that were significant in the mvMAPIT analysis.
#'
#' @format ## `mvmapit_data`
#' A nested list containing tibble data frames:
#' \describe{
#'   \item{mvmapit}{mvmapit return object; named list containing tibbles `pvalues`, `pves`, and `duration`.}
#'   \item{fisher}{Tibble containing fisher combiend p-values of the mvmapit data.}
#'   \item{exhaustive_search}{A dataframe containing the p-values of an exhaustive search together with the analysed interaction pair.}
#' }
#' @source data-raw/mvmapit_on_simulated_data.R
"mvmapit_data"
