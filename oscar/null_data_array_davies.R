library('MAPIT')
library('mvMAPIT')
library('foreach')


args = commandArgs(trailingOnly=TRUE)
datadir  <- Sys.getenv(c("SIMULATIONS_DIR"))
data_files <- list.files(path = file.path(datadir, "data/control"), pattern = 'null_H[0-9]_R[0-9]{1,2}_P10_.*rds')
mvmapit_results <- list()
task_id <- as.integer(args[1]) # Sys.getenv("SLURM_ARRAY_TASK_ID")
print(paste("Task ID", task_id))
f <- data_files[task_id]
print(paste("Data file", f))
st=format(Sys.time(), "%Y-%m-%d_%H:%M_")
mvmapit_out <- paste("mvMAPIT_davies_",st, f, sep = "")
print(mvmapit_out)
simulated <- readRDS(file.path(datadir, 'data/control', f))
mvmapit_results <- foreach(s=simulated[1:10]) %do% {
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
    hybrid = FALSE,
    threshold = 0.05,
    test = "davies",
    cores = 32,
    variantIndex = seq_len(ncol(X)),
    phenotypeCovariance = "identity",
    logLevel = "DEBUG"
  )
}
print('Save data.')
saveRDS(mvmapit_results, file=file.path(datadir, 'data/control/davies', mvmapit_out))

print('Finished mvMAPIT.')
