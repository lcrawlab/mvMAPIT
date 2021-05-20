library('mvMAPIT')
library('foreach')

args = commandArgs(trailingOnly=TRUE)
task_id <- as.integer(args[1]) # Sys.getenv("SLURM_ARRAY_TASK_ID")
print(paste("Task ID", task_id))
datadir  <- Sys.getenv(c("SIMULATIONS_DIR"))
data_files <- list.files(path = file.path(datadir, "data/epistasis"), pattern = 'H[0-9]_R[0-9]{1,2}_.*rds')
mvmapit_results <- list()
f <- data_files[task_id]
print(paste("Data file", f))
st=format(Sys.time(), "%Y-%m-%d_%H:%M_")
mvmapit_out <- paste("mvMAPIT_",st, f, sep = "")
print(mvmapit_out)
simulated <- readRDS(file.path(datadir, 'data/epistasis', f))
mvmapit_results <- foreach(s=simulated) %do% {
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
    cores = 32,
    variantIndex = NULL,
    phenotypeCovariance = "identity",
    logLevel = "DEBUG"
  )
}
saveRDS(mvmapit_results, file=file.path(datadir, 'data/epistasis/out', mvmapit_out))

print('Finished mvMAPIT.')
