library(mvMAPIT) # install mvMAPIT prior to running this file
library(parallel)
library(stats)
library(dplyr)
library(RcppAlgos)
load('../data/simulated_data.rda')

X <- simulated_data$genotype
y <- simulated_data$trait


# run mvMAPIT
cores <- detectCores()

mvmapit_hybrid <- mvmapit(
  t(X),
  t(y),
  cores = cores,
  logLevel = "DEBUG"
)

fisher <- fishers_combined(mvmapit_hybrid$pvalues)

# exhaustive search for p-values
thresh <- 0.05 / nrow(X) # Set a significance threshold

significant_snps <-  fisher %>%
    filter(p < thresh) # Call only marginally significant SNPs
pairs <- NULL
if (nrow(significant_snps) > 1) {
  pairnames <- comboGeneral(significant_snps$id, 2) # Generate unique pairs of SNP names; for length(names) = n, the result is a (n * (n-1)) x 2 matrix with one row corresponding to a pair
  for (k in seq_len(nrow(pairnames))) {
    fit <- lm(y ~ X[, pairnames[k, 1]]:X[, pairnames[k, 2]])
    p_value1 <- coefficients(summary(fit))[[1]][2, 4]
    p_value2 <- coefficients(summary(fit))[[2]][2, 4]
    tib <- tibble::tibble(
            x = p_value1,
            y = p_value2,
            u = pairnames[k, 1],
            v = pairnames[k, 2]
    )
    pairs <- bind_rows(pairs, tib)
  }
}

colnames(pairs) <- c(colnames(y), "var1", "var2")


# export data
mvmapit_data <- list(
  mvmapit = mvmapit_hybrid,
  fisher = fisher,
  exhaustive_search = pairs
)

usethis::use_data(mvmapit_data, overwrite = TRUE)
