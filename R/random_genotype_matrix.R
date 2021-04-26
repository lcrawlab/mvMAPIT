#' Random Genotype Matrix.
#'
#' A dataset containing 2938 samples with 5747 SNPs.
#'
#' @format A data frame with 2938 rows and 5747 variables:
#' \describe{
#'   \item{sample_name}{Unique name of each simulated sample.}
#'   \item{SNPs}{Matrix entries take values 0,1,2. 0 for homozygous in major allele, 1 for heterozygous, 2 for homozygous in minor allele.}
#' }
#' @source data-raw/random_data.R
"random_genotype_matrix"
