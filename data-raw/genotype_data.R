set.seed(1234)
library(dplyr)

n_samples <- 2938
n_snp <- 5747

sample_names <- seq_len(n_samples) %>% sprintf(fmt = "id%04d")
snp_names <- seq_len(n_snp) %>% sprintf(fmt = "snp%04d")

random_genotype_vec <- sample(0:2, n_samples * n_snp, replace = TRUE)
genotype_data <- matrix(random_genotype_vec,
  nrow = n_samples,
  ncol = n_snp,
)

colnames(genotype_data) <- snp_names
rownames(genotype_data) <- sample_names

usethis::use_data(genotype_data, overwrite = TRUE)
