library(tidyr)
library(tidyverse)
library(dplyr)
library(mvMAPIT)
library(ggplot2)
library(qqplotr)

# for alternative data -- need to toggle to find a way to run for null data or could write separate script

mvmapit_simulation <- function(data_type = c("null", "alternative")) {
    
    df <- NULL
    df_statistics <- data.frame(trial = c(), correlation_p = c()) #add heretibility and hypo 
    
    if (data_type == "null"){
        n_pleiotropic = 0
        n_trait_specific = 0
        print("null")
    } else {
        n_pleiotropic = 10
        n_trait_specific = 10
        print("alternative")
    }
    
    for (x in 1:5) {
        
        simulated_data <- simulate_traits(
            genotype_data,
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
            t(simulated_data$genotype[1:100, 1:100]),
            t(simulated_data$trait[1:100, ]),
            test = "normal", #could add into method inputs the type of test we want to run - matcharg
            skipProjection = FALSE
        )
        
        mvmapit_projection <- mvmapit(
            t(simulated_data$genotype[1:100, 1:100]),
            t(simulated_data$trait[1:100, ]),
            test = "normal", 
            skipProjection = TRUE
        )
        
        #store data in tidyr where p-value lists are diff for each trial
        df <- rbind(df, mvmapit$pvalues %>% mutate(trial = x, projections = "False")) 
        df <- rbind(df, mvmapit_projection$pvalues %>% mutate(trial = x, projections = "True")) 
        
        #print(cor(mvmapit$pvalues$p, mvmapit_projection$pvalues$p))
        df_statistics <- rbind(df_statistics, c(trial = x, correlation_p = cor(mvmapit$pvalues$p, mvmapit_projection$pvalues$p)))
    }
  
    #compare p-values af entire 
    

    #expectation under null -- pvalue uniform distribution (compare between proj and not)
    dp <- list(rate=log(10))
    di <- "exp"
    de <- FALSE # enabling the detrend option
    
    # TODO: change causal snp facet labels
    gg <- df %>% ggplot(mapping = aes(
        sample = -log10(p)
    )) +
        theme_bw() +
        stat_qq_band(distribution = di,
                     dparams = dp,
                     detrend = de,
                     alpha = 0.5) +
        stat_qq_line(distribution = di, dparams = dp, detrend = de) +
        stat_qq_point(distribution = di, dparams = dp, detrend = de) +
        theme(legend.position = "none") +
        labs(x = bquote("Theoretical Quantiles " -log[10](p)),
             y = bquote("Sample Quantiles " -log[10](p))) + 
        facet_wrap(~projections) # check this one
    show(gg)
    
    #evaluate distribution of p-value correlations
    colnames(df_statistics) <- c('trial', 'correlation_p')
    
    box<- ggplot(df_statistics, aes(y = correlation_p)) +
        geom_boxplot() + scale_x_discrete() +
        labs(title = "Projection correlation distribution",
             y = "correlation coefficent")
    show(box)
    return(box)
}   

box <- mvmapit_simulation("null")
show(box)




