library('MAPIT')
library('mvMAPIT')
library('foreach')

datadir  <- Sys.getenv(c("SIMULATIONS_DIR"))
data_files <- list.files(path = file.path(datadir, "data/control"), pattern = 'null_H[0-9]_R[0-9]{1,2}_.*rds')
mvmapit_results <- list()
for(f in data_files[1]) {
  st=format(Sys.time(), "%Y-%m-%d_%H:%M")
  mvmapit_out <- paste("mvMAPIT_",st, ".rds", sep = "")
  mapit_out <- paste("MAPIT_",st, ".rds", sep = "")
  simulated <- readRDS(file.path(datadir, 'data/control', f))
  mapit_results <- foreach(s=simulated[1]) %do% {
    y <- s$phenotype[,1]
    X <- s$genotype
    maf <- colMeans(X) / 2
    X <- X[, (maf > 0.01)]
    MAPIT(
      t(X),
      y
    )
  }
  mvmapit_results <- foreach(s=simulated[1]) %do% {
    Y <- s$phenotype[,1]
    X <- s$genotype
    maf <- colMeans(X) / 2
    X <- X[, (maf > 0.01)]
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
      logLevel = "DEBUG"
    )
  }
  saveRDS(mvmapit_results, file=file.path(datadir, 'data/control/out', mvmapit_out))
  saveRDS(mapit_results, file=file.path(datadir, 'data/control/out', mapit_out))
}
print('Finished mvMAPIT.')
