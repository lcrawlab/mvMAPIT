#!/usr/bin/env Rscript
# Benchmark script for mvMAPIT performance testing
# Creates baseline and compares optimization improvements

library('mvMAPIT')

cat("=== mvMAPIT Performance Benchmark ===\n\n")

# Benchmark configuration
configs <- list(
  small = list(ind = 50, nsnp = 50, name = "Small (50×50)"),
  medium = list(ind = 100, nsnp = 100, name = "Medium (100×100)"),
  large = list(ind = 200, nsnp = 200, name = "Large (200×200)")
)

methods <- c("normal", "davies", "hybrid")

# Function to generate test data
generate_test_data <- function(ind, nsnp, seed = 12345) {
  set.seed(seed)
  H2 <- 0.6
  rho <- 0.5
  maf <- 0.05 + 0.45*runif(nsnp)
  X <- (runif(ind*nsnp) < maf) + (runif(ind*nsnp) < maf)
  X <- matrix(as.double(X), ind, nsnp, byrow = TRUE)

  sim <- simulate_traits(X,
                        n_causal = min(80, floor(nsnp * 0.8)),
                        H2 = H2,
                        rho = rho,
                        logLevel = 'WARN',
                        seed = seed)

  # Standardize X
  Xmean <- apply(X, 2, mean)
  Xsd <- apply(X, 2, sd)
  X <- t((t(X) - Xmean) / Xsd)

  list(X = X, Y = sim$trait, ind = ind, nsnp = nsnp)
}

# Function to benchmark a single configuration
benchmark_config <- function(config, test_method, cores = 4) {
  cat(sprintf("  Testing %s with test=%s...", config$name, test_method))

  data <- generate_test_data(config$ind, config$nsnp)

  # Warm-up run (not timed)
  invisible(gc())

  # Timed run
  start_time <- Sys.time()
  mem_before <- gc()[2, 2]  # Total MB before

  result <- mvmapit(t(data$X),
                   t(data$Y),
                   cores = cores,
                   test = test_method,
                   logLevel = 'ERROR')

  mem_after <- gc()[2, 2]  # Total MB after
  end_time <- Sys.time()

  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  time_per_variant <- elapsed / data$nsnp
  mem_used <- mem_after - mem_before

  cat(sprintf(" %.2fs (%.3fs/variant)\n", elapsed, time_per_variant))

  list(
    config = config$name,
    test_method = test_method,
    ind = data$ind,
    nsnp = data$nsnp,
    total_time = elapsed,
    time_per_variant = time_per_variant,
    mem_used_mb = mem_used,
    timestamp = Sys.time()
  )
}

# Run all benchmarks
run_all_benchmarks <- function(cores = 4) {
  results <- list()

  for (config_name in names(configs)) {
    config <- configs[[config_name]]
    cat(sprintf("\n%s:\n", config$name))

    for (test_method in methods) {
      result <- benchmark_config(config, test_method, cores)
      results[[length(results) + 1]] <- result
    }
  }

  # Convert to data frame
  df <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      config = x$config,
      test_method = x$test_method,
      ind = x$ind,
      nsnp = x$nsnp,
      total_time = x$total_time,
      time_per_variant = x$time_per_variant,
      mem_used_mb = x$mem_used_mb,
      timestamp = as.character(x$timestamp),
      stringsAsFactors = FALSE
    )
  }))

  df
}

# Detect available cores
available_cores <- parallel::detectCores()
cores_to_use <- min(4, available_cores)

cat(sprintf("Using %d cores (of %d available)\n", cores_to_use, available_cores))

# Run benchmarks
results <- run_all_benchmarks(cores = cores_to_use)

# Print summary table
cat("\n=== Benchmark Results ===\n\n")
print(results[, c("config", "test_method", "total_time", "time_per_variant", "mem_used_mb")])

# Save results
output_file <- "tests/benchmark-results.rds"
baseline_file <- "tests/benchmark-baseline.rds"

# If baseline doesn't exist, create it
if (!file.exists(baseline_file)) {
  saveRDS(results, baseline_file)
  cat(sprintf("\n✓ Baseline saved to %s\n", baseline_file))
} else {
  # Compare to baseline
  baseline <- readRDS(baseline_file)

  cat("\n=== Comparison to Baseline ===\n\n")

  for (i in 1:nrow(results)) {
    current <- results[i, ]
    matching_baseline <- baseline[baseline$config == current$config &
                                  baseline$test_method == current$test_method, ]

    if (nrow(matching_baseline) > 0) {
      speedup <- matching_baseline$total_time / current$total_time
      improvement <- (1 - (current$total_time / matching_baseline$total_time)) * 100

      cat(sprintf("%s (%s): %.2fs -> %.2fs (%.1f%% %s, %.2fx speedup)\n",
                 current$config,
                 current$test_method,
                 matching_baseline$total_time,
                 current$total_time,
                 abs(improvement),
                 if (improvement > 0) "faster" else "slower",
                 speedup))
    }
  }
}

# Always save current results
saveRDS(results, output_file)
cat(sprintf("\n✓ Current results saved to %s\n", output_file))

cat("\n=== Benchmark Complete ===\n")
