library(tidyr)
library(tidyverse)
library(dplyr)
library(mvMAPIT)
library(ggplot2)

# for alternative data -- need to toggle to find a way to run for null data or could write separate script

mvmapit_simulation <- function(data_type = c("null", "alternative")) {
    
    total <- NULL
    df_statistics <- data.frame(trial = c(), correlation_p = c()) #add heretibility and hypo 
    genotype_data <- readRDS("/oscar/data/lcrawfo1/sli347/biobank_1_10000.rds")
    
    if (data_type == "null"){
        n_pleiotropic = 0
        n_trait_specific = 0
    } else {
        n_pleiotropic = 10
        n_trait_specific = 10
    }
    
    for (x in 1:5) {

        sample_ids <- sample(1:10000, 5000, replace = FALSE)

        simulated_data <- simulate_traits(
            genotype_data[sample_ids, 1:5000],
            n_causal = 1000,
            n_trait_specific,
            n_pleiotropic,
            H2 = 0.6, #error = 1-H^2 -- could vary
            d = 1,
            rho = 0.2, #decrease number means easier to find snp
            marginal_correlation = 0.2,
            epistatic_correlation = 0.2,
            group_ratio_trait = 1,
            group_ratio_pleiotropic = 1,
            maf_threshold = 0.01,
            seed = sample(c(1:1e6), 1), #random generate
            logLevel = "INFO",
            logFile = NULL
        )

        #run on smaller data set -- to test
        #might want to randomize columns


        mvmapit <- mvmapit(
            t(simulated_data$genotype),
            t(simulated_data$trait),
            test = "normal", #could add into method inputs the type of test we want to run - matcharg
            skipProjection = FALSE
        )

        mvmapit_projection <- mvmapit(
            t(simulated_data$genotype),
            t(simulated_data$trait),
            test = "normal",
            skipProjection = TRUE
        )

        #store data in tidyr where p-value lists are diff for each trial
        total <- rbind(total, mvmapit$pvalues %>% mutate(trial = x, projection = "False"))
        total <- rbind(total, mvmapit_projection$pvalues %>% mutate(trial = x, projection = "True"))

        #print(cor(mvmapit$pvalues$p, mvmapit_projection$pvalues$p))
        df_statistics <- rbind(df_statistics, c(trial = x, correlation_p = cor(mvmapit$pvalues$p, mvmapit_projection$pvalues$p)))
    }

    #compare p-values af entire
    #expectation under null -- pvalue uniform distribution (compare between proj and not)
    ntest <- total %>%
        group_by(projection) %>%
        summarize(n = n())

    ci <- 0.95
    n_snps <- ntest$n[1]
    clower = -log10(qbeta(
        p = (1 - ci) / 2,
        shape1 = seq(n_snps),
        shape2 = rev(seq(n_snps))
    ))
    cupper = -log10(qbeta(
        p = (1 + ci) / 2,
        shape1 = seq(n_snps),
        shape2 = rev(seq(n_snps))
    ))

    df1 <- total %>%
        group_by(projection) %>%
        arrange(p, .by_group = TRUE) %>%
        mutate(
            cupper = cupper,
            clower = clower,
            expected= -log10(ppoints(n_snps)),
            observed = -log10(p)
        )

    gg <- df1 %>%
        ggplot(mapping= aes(x = expected, y = observed, color = projection)) +
        theme_bw() +
        facet_wrap(~projection) +
        geom_ribbon(aes(ymax = cupper, ymin = clower),
                    fill = "#999999", alpha = 0.5,linewidth =0.25
        ) +
        geom_point() +
        geom_segment(
            data = . %>% filter(expected == max(expected)),
            aes(x = 0, xend = expected, y = 0, yend = expected),
            linewidth = 1, alpha = 1, lineend = "round"
        ) +
        geom_hline(yintercept=-log10(0.05 / ntest$n[1]), linetype='dashed', col = 'black') +
        labs(x = bquote("Theoretical Quantiles " -log[10](p)),
             y = bquote("Sample Quantiles " -log[10](p))) +
        theme(legend.position = "none") +
        scale_color_manual(values = c("SteelBlue", "Indianred", "lightblue3", "#d79c9c", "darkgreen"))

    ggsave("/oscar/data/lcrawfo1/sli347/repeated_qqplot3.png", plot = gg, width = 6, height = 4, unit = "in", dpi = 300)


    #evaluate distribution of p-value correlations
    colnames(df_statistics) <- c('trial', 'correlation_p')

    box <- ggplot(df_statistics, aes(y = correlation_p)) +
        geom_boxplot() + scale_x_discrete() +
        labs(title = "Projection correlation distribution",
             y = "correlation coefficent")
    ggsave("/oscar/data/lcrawfo1/sli347/repeated_boxplot.png", plot = box, width = 6, height = 4, unit = "in", dpi = 300)

}   

mvmapit_simulation("null")



