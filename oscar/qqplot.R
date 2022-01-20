library('foreach')
library('gap')

print('Create qq-plots')
main_dir  <- Sys.getenv(c("DATA_DIR"))
data_dir <- file.path(main_dir, c("normal/control/combinatorial")) #, "davies/control"))
data_files <- foreach(d=data_dir, .combine = c) %do% {
  file_names <- list.files(path = d, pattern = 'mvMAPIT_.*rds', include.dirs = TRUE)
  foreach(n=file_names, .combine = c) %do% {
    file.path(d, n)
  }
}
for(f in data_files) {
  print(f)
  mapit <- readRDS(f)
  pvalues <- foreach(s=mapit, .combine = c) %do% {
    p <- as.numeric(s$pvalues)
    p <- p[!is.na(p)]
  }
  png(sub('\\.rds', '_qqplot.png', f), width=300, height=450)
  qqunif(pvalues, ci=TRUE)
  dev.off()
}
print('Finished qq-plots.')
