set.seed(1234)
library(dplyr)

n_samples <- 2938
n_snp <- 5747

sample_names <- seq_len(n_samples) %>% sprintf(fmt = "id%04d")
snp_names <- seq_len(n_snp) %>% sprintf(fmt = "snp%04d")

random_genotype_data <- sample(0:2, n_samples * n_snp, replace = TRUE)
random_genotype_matrix <- matrix(random_genotype_data,
  nrow = n_samples,
  ncol = n_snp,
)

colnames(random_genotype_matrix) <- snp_names
rownames(random_genotype_matrix) <- sample_names

usethis::use_data(random_genotype_matrix, overwrite = TRUE)
