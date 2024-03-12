# Setup a trial specification using a binary, binomially
# distributed, undesirable outcome
binom_trial <- setup_trial_binom(
    arms = c("Arm A", "Arm B", "Arm C"), 
    true_ys = c(0.25, 0.20, 0.30),
    # Minimum allocation of 15% in all arms
    min_probs = rep(0.15, 3),
    data_looks = seq(from = 300, to = 2000, by = 100),
    # Stop for equivalence if > 90% probability of
    # absolute differences < 5 percentage points
    equivalence_prob = 0.9,
    equivalence_diff = 0.05,
    soften_power = 0.5 # Limit extreme allocation ratios
)


# Run 10 simulations with a specified random base seed
res <- run_trials(binom_trial, n_rep = 10, base_seed = 12345) #how to implement simulation and ohter code?? 
# in the mix?

# See ?extract_results, ?check_performance, ?summary and ?print for details
# on extracting resutls, summarising and printing

p-values <- list()

for (x in 1:100) {
    file_name <- alt_data + str(x)
    alt_data <- simulate_traits(
        genotype_data,
        n_causal = 1000,
        n_trait_specific = runif(n=1, min=1, max=20), #random 
        n_pleiotropic = 0, 
        H2 = 0.6,
        d = 1,
        rho = runif(1), #random
        marginal_correlation = 0.2, 
        epistatic_correlation = 0.2,
        group_ratio_trait = 1,
        group_ratio_pleiotropic = 1,
        maf_threshold = 0.01,
        seed = 67132,
        logLevel = "INFO",
        logFile = NULL
    )
    
    # Save an object to a file
    saveRDS(alt_data, file = file_name + ".rds")
    # Restore the object
    my_data <- readRDS(file = "~/mvMAPIT/" + file_name + ".rds")
    
    mvmapit_normal_alt <- mvmapit(
        t(null_data$genotype),
        t(null_data$trait),
        test = "normal"
    )
    fisher <- fishers_combined(mvmapit_normal_alt$pvalues) #do we need this part
    
    
}
