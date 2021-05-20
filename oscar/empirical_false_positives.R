library('foreach')

simulations_dir  <- Sys.getenv(c("SIMULATIONS_DIR"))
data_dir <- "data/control/davies"
results_dir <- "data/control/davies/results"
data_files <- list.files(path = file.path(simulations_dir, data_dir), pattern = '.*rds')
thresholds <- c(0.05, 0.01, 0.001)
results_mvmapit <- matrix(nrow = 0, ncol = 16)
for(f in data_files) {
  mapit <- readRDS(file.path(simulations_dir, data_dir, f))
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
colnames(df) <- c('N', 'H2', 'rho', 'P', 't05', 't01', 't001', 'dt05', 'dt01', 'dt001', 'K', 'M', 'kron', 'q', 'S', 'p')
print(results_mvmapit)
write.csv(df, file=file.path(simulations_dir, results_dir, 'mvMAPIT.csv'), row.names=FALSE)
saveRDS(results_mvmapit, file=file.path(simulations_dir, results_dir, 'mvMAPIT.rds'))
print('Finished mvMAPIT.')
