library('mvMAPIT')
library('foreach')

datadir  <- Sys.getenv(c("MVMAPIT_DATA_DIR"))
data_files <- list.files(path = file.path(datadir, "data/control"), pattern = 'null_H[0-9]_R[0-9]{1,2}_.*rds')
results <- list()
for(f in data_files[1]) {
  simulated <- readRDS(file.path(datadir, 'data/control', f))
  results[[f]] <- foreach(s=simulated[1]) %do% {
    #ind <- sample(seq_len(nrow(s$genotype)), 500, replace = FALSE)
    Y <- s$phenotype[, 1:2]
    X <- s$genotype
    #X <- X[, -as.numeric(which(apply(X, 2, var) == 0))]
    MvMAPIT(
      t(X),
      t(Y),
      W = NULL,
      C = NULL,
      hybrid = TRUE,
      threshold = 0.05,
      test = "normal",
      cores = 1,
      variantIndex = NULL,
      phenotypeCovariance = "identity",
      logLevel = "DEBUG",
      #logFile = 'eval_type_1.log'
    )
  }
}
saveRDS(results, file=file.path(datadir, 'data/control/out', 'results.rds'))
print('Finished mvMAPIT.')
