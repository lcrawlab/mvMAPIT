#' Fisher's combine method p-value computation
#'
#' This function takes in p-values to combine them via the Fisher's method.
#'
#' @param pvalues Vector with p-values to combine
#' @return Scalar Fisher's combined p-value
#' @noRd
#' @importFrom stats pchisq
sumlog <- function(pvalues) {
    df <- 2 * length(pvalues)
    fisherp <- pchisq(
        -2 * sum(log(pvalues)),
        df, lower.tail = FALSE
    )
    return(fisherp)
}

#' Fisher's combine method on mvmapit return
#'
#' This function takes in the p-values tibble that mvmapit returned. It then
#' computes the combined p-values grouped by variant id.
#'
#' @param pvalues Tibble with p-values from mvmapit function call.
#' @param group_col String that denotes column by which to group and combine
#' p-values.
#' @param p_col String that denotes p-value column.
#' @return A Tibble with the combined p-values.
#' @examples
#' set.seed(837)
#' p <- 200
#' n <- 100
#' d <- 2
#' X <- matrix(
#'     runif(p * n),
#'     ncol = p
#' )
#' Y <- matrix(
#'     runif(d * n),
#'     ncol = d
#' )
#' mapit <- mvmapit(
#'     t(X),
#'     t(Y),
#'     test = "normal", cores = 1, logLevel = "INFO"
#' )
#' fisher <- fishers_combined(mapit$pvalues)
#' @export
#' @import dplyr
fishers_combined <- function(pvalues, group_col = "id", p_col = "p") {
    pvalues %>%
        group_by(.data[[group_col]]) %>%
        summarize(p = sumlog(.data[[p_col]])) %>%
        mutate(trait = "fisher") %>%
        relocate(all_of(group_col), all_of("trait"), all_of("p"))
}


#' Harmonic mean p combine method on mvmapit return
#'
#' This function takes in the p-values tibble that mvmapit returned. It then
#' computes the combined p-values grouped by variant id.
#'
#' @param pvalues Tibble with p-values from mvmapit function call. Grouping is
#' based on the column named "id"
#' @param group_col String that denotes column by which to group and combine
#' p-values.
#' @param p_col String that denotes p-value column.
#' @return A Tibble with the combined p-values.
#' @examples
#' set.seed(837)
#' p <- 200
#' n <- 100
#' d <- 2
#' X <- matrix(
#'     runif(p * n),
#'     ncol = p
#' )
#' Y <- matrix(
#'     runif(d * n),
#'     ncol = d
#' )
#' mapit <- mvmapit(
#'     t(X),
#'     t(Y),
#'     test = "normal", cores = 1, logLevel = "INFO"
#' )
#' harmonic <- harmonic_combined(mapit$pvalues)
#' @export
#' @import harmonicmeanp
#' @import dplyr
harmonic_combined <- function(pvalues, group_col = "id", p_col = "p") {
    pvalues %>%
        group_by(.data[[group_col]]) %>%
        summarize(p = as.numeric(hmp.stat(.data[[p_col]]))) %>%
        mutate(trait = "harmonic") %>%
        relocate(all_of(group_col), all_of("trait"), all_of("p"))
}


#' Cauchy p combine method on mvmapit return
#'
#' This function takes in the p-values tibble that mvmapit returned. It then
#' computes the combined p-values grouped by variant id.
#'
#' @param pvalues Tibble with p-values from mvmapit function call. Grouping is
#' based on the column named "id"
#' @param group_col String that denotes column by which to group and combine
#' p-values.
#' @param p_col String that denotes p-value column.
#' @return A Tibble with the combined p-values.
#' @examples
#' set.seed(837)
#' p <- 200
#' n <- 100
#' d <- 2
#' X <- matrix(
#'     runif(p * n),
#'     ncol = p
#' )
#' Y <- matrix(
#'     runif(d * n),
#'     ncol = d
#' )
#' mapit <- mvmapit(
#'     t(X),
#'     t(Y),
#'     test = "normal", cores = 1, logLevel = "INFO"
#' )
#' cauchy <- cauchy_combined(mapit$pvalues)
#' @export
#' @import dplyr
cauchy_combined <- function(pvalues, group_col = "id", p_col = "p") {
    pvalues %>%
        group_by(.data[[group_col]]) %>%
        summarize(p = as.numeric(cauchy_p(.data[[p_col]]))) %>%
        mutate(trait = "cauchy") %>%
        relocate(all_of(group_col), all_of("trait"), all_of("p"))
}


#' An analytical p-value combination method using the Cauchy distribution
#'
#' The \code{cauchy_p} function takes in a numeric vector of p-values, a numeric
#' vector of non-negative weights, and return the aggregated p-value using Cauchy method.
#' @param pvals a numeric vector of p-values, where each of the element is
#' between 0 to 1, to be combined.
#' @param weights a numeric vector of non-negative weights. If \code{NULL}, the
#' equal weights are assumed (default = NULL).
#' @return the aggregated p-value combining p-values from the vector \code{pvals}.
#' @examples
#' pvalues <- c(2e-02, 4e-04, 0.2, 0.1, 0.8)
#' cauchy_p(pvals = pvalues)
#' @references Liu, Y., & Xie, J. (2020). Cauchy combination test: a powerful test
#' with analytic p-value calculation under arbitrary dependency structures.
#' \emph{Journal of the American Statistical Association}, \emph{115}(529), 393-402.
#' (\href{https://doi.org/10.1080/01621459.2018.1554485}{pub})
#' @references Liu, Y., et al. (2019). Acat: A fast and powerful p value combination
#' method for rare-variant analysis in sequencing studies.
#' \emph{The American Journal of Human Genetics}, \emph{104}(3), 410-421.
#' (\href{https://doi.org/10.1016/j.ajhg.2019.01.002}{pub})
#' @references Li, X., et al. (2020). Dynamic incorporation of multiple in silico
#' functional annotations empowers rare variant association analysis of large
#' whole-genome sequencing studies at scale. Nature genetics, 52(9), 969-983.
#' @references https://github.com/xihaoli/STAAR
#' @noRd
#' @importFrom stats pcauchy
cauchy_p <- function(pvals, weights=NULL){
  #### check if there is NA
  if(sum(is.na(pvals)) > 0){
    stop("Cannot have NAs in the p-values!")
  }

  #### check if all p-values are between 0 and 1
  if((sum(pvals<0) + sum(pvals>1)) > 0){
    stop("All p-values must be between 0 and 1!")
  }

  #### check if there are p-values that are either exactly 0 or 1.
  is.zero <- (sum(pvals==0)>=1)
  is.one <- (sum(pvals==1)>=1)
  if(is.zero && is.one){
    stop("Cannot have both 0 and 1 p-values!")
  }
  if(is.zero){
    return(0)
  }
  if(is.one){
    warning("There are p-values that are exactly 1!")
    return(1)
  }

  #### check the validity of weights (default: equal weights) and standardize them.
  if(is.null(weights)){
    weights <- rep(1/length(pvals),length(pvals))
  }else if(length(weights)!=length(pvals)){
    stop("The length of weights should be the same as that of the p-values!")
  }else if(sum(weights < 0) > 0){
    stop("All the weights must be positive!")
  }else{
    weights <- weights/sum(weights)
  }

  #### check if there are very small non-zero p-values
  is.small <- (pvals < 1e-16)
  if (sum(is.small) == 0){
    cct.stat <- sum(weights*tan((0.5-pvals)*pi))
  }else{
    cct.stat <- sum((weights[is.small]/pvals[is.small])/pi)
    cct.stat <- cct.stat + sum(weights[!is.small]*tan((0.5-pvals[!is.small])*pi))
  }

  #### check if the test statistic is very large.
  if(cct.stat > 1e+15){
    pval <- (1/cct.stat)/pi
  }else{
    pval <- 1-pcauchy(cct.stat)
  }
  return(pval)
}
