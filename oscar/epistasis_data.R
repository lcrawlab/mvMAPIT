library('MAPIT')
library('mvMAPIT')
library('foreach')

datadir  <- Sys.getenv(c("SIMULATIONS_DIR"))
data_files <- list.files(path = file.path(datadir, "data/epistasis"), pattern = 'H[0-9]_R[0-9]{1,2}_.*rds')
mvmapit_results <- list()
for(f in data_files) {
  st=format(Sys.time(), "%Y-%m-%d_%H:%M_")
  mvmapit_out <- paste("mvMAPIT_",st, f, sep = "")
  mapit_out <- paste("MAPIT_",st, f, sep = "")
  print(mvmapit_out)
  simulated <- readRDS(file.path(datadir, 'data/epistasis', f))
  mapit_results <- foreach(s=simulated[1]) %do% {
    y <- s$phenotype[,1]
    X <- s$genotype
    print(dim(X))
    X <- X[, which(apply(X, 2, var) != 0)]
    print(dim(X))
    maf <- colMeans(X) / 2
    X <- X[, (maf > 0.01)]
    MAPIT(
      t(X),
      y
    )
  }
  mvmapit_results <- foreach(s=simulated[1:3]) %do% {
    Y <- s$phenotype
    X <- s$genotype
    print(dim(X))
    X <- X[, which(apply(X, 2, var) != 0)]
    print(dim(X))
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
  saveRDS(mvmapit_results, file=file.path(datadir, 'data/epistasis/out', mvmapit_out))
  saveRDS(mapit_results, file=file.path(datadir, 'data/epistasis/out', mapit_out))
}
print('Finished mvMAPIT.')
