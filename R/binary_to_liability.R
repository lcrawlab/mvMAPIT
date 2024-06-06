#' Convert binary traits to liabilities
#'
#' This function implements the conversion of binary traits to liabilties as
#' proposed in the LT-MAPIT model (Crawford and Zhou 2018,
#' https://doi.org/10.1101/374983).
#' To run LT-MAPIT (MAPIT on case-control traits), convert the binary traits to
#' liabilities using this function and pass the liabilities to mvmapit as trait.
#'
#' @param case_control_trait Case-control trait encoded as binary trait with 0 as control or 1 as case.
#' @param prevalence Case prevalence between 0 and 1. Proportion of cases in the population.
#' @returns A trait vector of same length as y with case-control indicators converted
#' to liabilties.
#'
#' @export
#' @importFrom truncnorm rtruncnorm
#' @importFrom stats qnorm
#' @import checkmate
binary_to_liability <- function(case_control_trait, prevalence) {
  coll <- makeAssertCollection()
  assertInteger(
    case_control_trait,
    lower = 0,
    upper = 1,
    add = coll
  )
  assertDouble(prevalence,
               lower = 0,
               upper = 1,
               add = coll)
  reportAssertions(coll)
  
  n_cases = sum(case_control_trait == 1, na.rm = T)
  n_controls = sum(case_control_trait == 0, na.rm = T)
  threshold = qnorm(1 - prevalence, mean = 0, sd = 1)
  
  liabilities = rep(NA, length(case_control_trait))
  liabilities[!is.na(case_control_trait) &
                case_control_trait == 0] = rtruncnorm(n_controls, b = threshold)
  liabilities[!is.na(case_control_trait) &
                case_control_trait == 1] = rtruncnorm(n_cases, a = threshold)
  
  return(liabilities)
}
