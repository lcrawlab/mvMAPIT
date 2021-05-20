library('foreach')
library('gap')

simulations_dir  <- Sys.getenv(c("DATA_DIR"))
args = commandArgs(trailingOnly=TRUE)
data_dir <- "normal/control"
data_files <- list.files(path = file.path(simulations_dir, data_dir), pattern = 'mvMAPIT.*rds')
for(f in data_files) {
  mapit <- readRDS(file.path(simulations_dir, data_dir, f))
  pvalues <- foreach(s=mapit, .combine = c) %do% {
    p <- as.numeric(s$pvalues)
    p <- p[!is.na(p)]
  }
  png(paste0(file.path(simulations_dir, data_dir), sub('\\.rds', '', f), '_qqplot.png'))
  qqunif(pvalues, ci=TRUE)
  dev.off()
}
print('Finished qq-plots.')
