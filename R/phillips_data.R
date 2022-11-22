#' Multivariate MAPIT analysis of binding affinities in broadly neutralizing antibodies.
#'
#' This data set contains the return object from the multivariate MAPIT method,
#' applied to two binding affinity traits for two broadly neutralizing antibodies.
#' It also contains the regression coefficients on the same data as reported by Phillips et al. (2021).
#' 
#' The antibody CR9114 was analyzed with influenza H1 and H3.
#' The antibody CR6261 was analyzed with influenza H1 and H9.
#' 
#' In the data, the p-values are computed for the test whether a given residue
#' position has a marginal epistatic effect on the binding affinities.
#' 
#' Phillips et al. (2021) Binding affinity landscapes constrain the evolution of broadly neutralizing anti-influenza antibodies. eLife 10:e71393
#'
#' @format
#' A named list containing tibble data frames:
#' \describe{
#'   \item{fisher}{Tibble containing among other columns the residue id, p-values, antibody species, trait of the Phillips data. Combined p-value with Fisher's method.}
#'   \item{harmonic}{Tibble containing among other columns the residue id, p-values, antibody species, trait of the Phillips data. Combined p-value with harmonic mean p method.}
#'   \item{regression}{Named list containing two tibbles containing regression coefficients as reported by Phillips et al.}
#' }
#' @source vignette/study-phillips-bnabs.Rmd
"phillips_data"
