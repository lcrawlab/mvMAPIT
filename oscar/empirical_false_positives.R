library('MAPIT')
library('mvMAPIT')
library('foreach')

simulations_dir  <- Sys.getenv(c("SIMULATIONS_DIR"))
data_dir <- "data/control/out"
results_dir <- "data/control/results"
data_files <- list.files(path = file.path(simulations_dir, data_dir), pattern = '.*rds')
thresholds <- c(0.05, 0.01, 0.001)
results_mvmapit <- matrix(nrow = 0, ncol = 7)
results_mapit <- matrix(nrow = 0, ncol = 7)
for(f in data_files) {
  mapit <- readRDS(file.path(simulations_dir, data_dir, f))
  pvalues <- foreach(s=mapit, .combine = c) %do% {
    as.numeric(s$pvalues)
  }
  type1 <- foreach(t=thresholds) %do% {
    sum(pvalues < t) / length(pvalues)
  }
  H <- as.numeric(gsub(".*H([0-9]+).*$", "\\1", f)) / 10
  R <- as.numeric(gsub(".*R([0-9]+).*$", "\\1", f)) / 10
  N <- as.numeric(gsub(".*N([0-9]+).*$", "\\1", f))
  P <- as.numeric(gsub(".*P([0-9]+).*$", "\\1", f))
  if (grepl('mvMAPIT', f, fixed = TRUE)) {
    results_mvmapit <- rbind(results_mvmapit, c(N, H, R, P, type1))
  } else {
    results_mapit <- rbind(results_mapit, c(N, H, R, P, type1))
  }
}
saveRDS(results_mvmapit, file=file.path(simulations_dir, results_dir, 'mvMAPIT.rds'))
saveRDS(results_mapit, file=file.path(simulations_dir, results_dir, 'MAPIT.rds'))
print('Finished mvMAPIT.')
