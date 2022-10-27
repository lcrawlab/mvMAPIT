#' Fisher's combine method p-value computation
#'
#' This function takes in p-values to combine them via the Fisher's method.
#'
#' @param pvalues Vector with p-values to combine
#' @return Scalar Fisher's combined p-value
#' @noRd
sumlog <- function(pvalues) {
    df <- 2 * length(pvalues)
    fisherp <- pchisq(
        -2 * sum(log(pvalues)),
        df, lower.tail = FALSE
    )
    return(fisherp)
}

#' Fisher's combine method on MvMAPIT return
#'
#' This function takes in the p-values tibble that MvMAPIT returned. It then
#' computes the combined p-values grouped by variant id.
#'
#' @param pvalues Tibble with p-values from MvMAPIT function call. Grouping is
#' based on the column named "id"
#' @return A Tibble with the combined p-values.
#' @export
fishers_combined <- function(pvalues) {
    pvalues %>%
        group_by(id) %>%
        summarize(p = sumlog(p)) %>%
        mutate(trait = "fisher") %>%
        relocate(id, trait, p)
}


#' Harmonic mean p combine method on MvMAPIT return
#'
#' This function takes in the p-values tibble that MvMAPIT returned. It then
#' computes the combined p-values grouped by variant id.
#'
#' @param pvalues Tibble with p-values from MvMAPIT function call. Grouping is
#' based on the column named "id"
#' @return A Tibble with the combined p-values.
#' @export
#' @import harmonicmeanp
harmonic_combined <- function(pvalues) {
    pvalues %>%
        group_by(id) %>%
        summarize(p = as.numeric(hmp.stat(p))) %>%
        mutate(trait = "harmonic") %>%
        relocate(id, trait, p)
}
