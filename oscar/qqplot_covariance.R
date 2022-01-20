library("foreach")
library("gap")

print("Create qq-plots for combinatorial mvMAPIT run")
main_dir  <- Sys.getenv(c("DATA_DIR"))
data_dir <- file.path(main_dir, c("normal/control/combinatorial")) #, "davies/control"))
data_files <- foreach(d = data_dir, .combine = c) %do% {
  file_names <- list.files(path = d, pattern = "mvMAPIT_.*rds", include.dirs = TRUE)
  foreach(n = file_names, .combine = c) %do% {
    file.path(d, n)
  }
}
for (f in data_files) {
  print(f)
  mapit <- readRDS(f)
  pvalues <- foreach(s = mapit, .combine = rbind) %do% {
    p <- s$pvalues
    row.has.na <- apply(p, 1, function(x){any(is.na(x))})
    p <- p[!row.has.na, ]
  }
  for (i in seq_len(ncol(pvalues))) {
    plot.name <- sub("\\.rds",
            paste0("_", sub("\\*", "", colnames(pvalues)[i]), "_qqplot.png"), f)
    png(plot.name, width = 300, height = 450)
    qqunif(pvalues[, i], ci = TRUE)
    dev.off()
    print(plot.name)
  }
}
print("Finished qq-plots.")
