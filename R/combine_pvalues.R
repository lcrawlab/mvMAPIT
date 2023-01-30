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
