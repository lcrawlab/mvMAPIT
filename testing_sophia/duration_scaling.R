library(tidyr)
library(tidyverse)
library(dplyr)
library(mvMAPIT)
library(ggplot2)

set.seed(1234)

n_samples <- c(1000, 2000, 5000, 10000)
n_snps <- c(100, 1000)

pairs <- crossing(n_samples, n_snps)
df <- NULL
df_concise <- NULL

for(i in 1:nrow(pairs)) {

    n_snp <- pairs$n_snps[i]
    n_samples <- pairs$n_samples[i]
    maf <- 0.05 + 0.45 * runif(n_snp)
    random_minor_allele_counts   <- (runif(n_samples * n_snp) < maf) + (runif(n_samples * n_snp) < maf)
    genotype_data <- matrix(random_minor_allele_counts,
                        nrow = n_samples,
                        ncol = n_snp,
                        byrow = TRUE,
    )

    sample_names <- seq_len(n_samples) %>% sprintf(fmt = "id%04d")
    snp_names <- seq_len(n_snp) %>% sprintf(fmt = "snp%04d")

    colnames(genotype_data) <- snp_names
    rownames(genotype_data) <- sample_names
    
    trait <- rnorm(n_samples)
    
    mvmapit_normal <- mvmapit(
        t(genotype_data),
        t(trait),
        test = "normal", #could add into method inputs the type of test we want to run - matcharg
        skipProjection = FALSE
    )
    
    mvmapit_normal_projection <- mvmapit(
        t(genotype_data),
        t(trait),
        test = "normal", 
        skipProjection = TRUE
    )
    
    df <- rbind(df, mvmapit_normal$duration %>% mutate(n_samples = n_samples, n_snp = n_snp, projection = "False")) 
    df <- rbind(df, mvmapit_normal_projection$duration %>% mutate(n_samples = n_samples, n_snp = n_snp, projection = "True")) 
    
    #df <- rbind(df, c(mvmapit_normal$duration$duration_ms, n_samples = n_samples,
                      #n_snps = n_snps, projections = "False"))
    
    #df <- rbind(df, c(mvmapit_normal_projection$duration$duration_ms, n_samples = n_samples, 
                      #n_snps = n_snps, projections = "True"))
    
}

#colnames(df) <- c("compute_covariances", "project_matrices", "vectorize", "compute_q", 
 #                         "compute_S", "compute_var_delta", "n_samples", "n_snps", "projections")

#df$total_duration <- rowSums(df[ , c(1, 2, 3, 4, 5, 6)], na.rm=TRUE)

result <- df %>% group_by(n_samples, n_snp, projection) %>% summarise(total = sum(duration_ms))

saveRDS(df, "/oscar/data/lcrawfo1/sli347/duration_scaling_df2.RDS")
saveRDS(result, "/oscar/data/lcrawfo1/sli347/duration_scaling2.RDS")
