library(doParallel)
library(Rcpp)
library(RcppAlgos)
library(RcppArmadillo)
library(RcppParallel)
library(CompQuadForm)

load('data/simulated_epistasis_data.rda')

X <- simulated_epistasis_data$genotype
y <- simulated_epistasis_data$phenotype

# run MAPIT
cores <- detectCores()

ptm <- proc.time() # Start clock
mapit_hybrid <- MAPIT(t(X), y, cores = cores)
proc.time() - ptm # Stop clock

hybrid.pvals <- mapit_hybrid$pvalues
names(hybrid.pvals) <- colnames(X)

ptm <- proc.time() # Start clock
mapit_normal <- MAPIT(
  t(X),
  y,
  hybrid = FALSE,
  test = "normal",
  cores = cores
)
proc.time() - ptm # Stop clock

normal.pvals1 <- mapit_normal$pvalues
names(normal.pvals1) <- colnames(X)

ptm <- proc.time() # Start clock
mapit_davies <- MAPIT(
  t(X),
  y,
  hybrid = FALSE,
  test = "davies",
  cores = cores
)
proc.time() - ptm # Stop clock

davies.pvals <- mapit_davies$pvalues
names(davies.pvals) <- colnames(X)

# exhaustive search for p-values
thresh <- 0.05 / length(hybrid.pvals) # Set a significance threshold
significant_snps <- hybrid.pvals[hybrid.pvals <= thresh] # Call only marginally significant SNPs
significant_snps <- significant_snps[!is.na(significant_snps)]
pairs <- c()
if (length(significant_snps) > 1) {
  pairnames <- comboGeneral(names(significant_snps), 2) # Generate unique pairs of SNP names; for length(names) = n, the result is a (n * (n-1)) x 2 matrix with one row corresponding to a pair
  for (k in seq_len(nrow(pairnames))) {
    fit <- lm(y ~ X[, pairnames[k, 1]]:X[, pairnames[k, 2]])
    p_value <- coefficients(summary(fit))[8]
    names(p_value) <- paste(pairnames[k, 1], pairnames[k, 2], sep = "-")
    pairs <- c(pairs, p_value)
  }
}