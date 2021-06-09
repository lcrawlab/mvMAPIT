library('foreach')

print('Empirical type 1 error computation')
main_dir  <- Sys.getenv(c("DATA_DIR"))
data_dir <- file.path(main_dir, c("normal/control", "davies/control"))
thresholds <- c(0.05, 0.01, 0.001)
for (d in data_dir) {
  results_mvmapit <- matrix(nrow = 0, ncol = 16)
  data_files <- list.files(path = d, pattern = 'mvMAPIT_.*rds')
  for(f in data_files) {
    print(f)
    mapit <- readRDS(file.path(d, f))
    pvalues <- foreach(s=mapit, .combine = rbind) %do% {
      p <- as.numeric(s$pvalues)
      p <- p[!is.na(p)]
      type1 <- foreach(t=thresholds, .combine = c) %do% {
        sum(p < t) / length(p)
      }
    }
    timings <- foreach(s=mapit, .combine = rbind) %do% {
      as.numeric(s$timings)
    }
    type1_mean <- apply(pvalues, 2, mean)
    type1_var <- apply(pvalues, 2, sd)
    times <- apply(timings, 2, mean)
    print(type1_mean)
    print(type1_var)
    H <- as.numeric(gsub(".*H([0-9]+).*$", "\\1", f)) / 10
    R <- as.numeric(gsub(".*R([0-9]+).*$", "\\1", f)) / 10
    N <- as.numeric(gsub(".*N([0-9]+).*$", "\\1", f))
    P <- as.numeric(gsub(".*P([0-9]+).*$", "\\1", f))
    results_mvmapit <- rbind(results_mvmapit, c(N, H, R, P, type1_mean, type1_var, times))
  }
  df <- as.data.frame(results_mvmapit)
  colnames(df) <- c('N', 'H2', 'rho', 'P', 'p05', 'p01', 'p001', 'dp05', 'dp01', 'dp001', 'K', 'M', 'kron', 'q', 'S', 'p')
  print(results_mvmapit)
  write.csv(df, file=file.path(d, 'results.csv'), row.names=FALSE)
  saveRDS(df, file=file.path(d, 'results.rds'))
}
print('Finished mvMAPIT.')
