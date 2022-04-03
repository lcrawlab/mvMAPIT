# Runs the Davies portion of the hypothesis testing
mvmapit_pvalues <- function(cpp_structure, X, accuracy) {
  num_combinations <- dim(cpp_structure$Eigenvalues)[2]
  p <- nrow(X)
  hybrid.pvals <- matrix(1, p, num_combinations)
  for (c in seq_len(num_combinations)) {
    hybrid.pvals[, c] <- hybrid_pvalues(cpp_structure$Est[, c],
                                        cpp_structure$Eigenvalues[, c,],
                                        accuracy)
  }
  names(hybrid.pvals) <- rownames(X)
  return(hybrid.pvals)
}

hybrid_pvalues <- function(point_estimates, eigenvalues, acc) {
  ### Apply Davies Exact Method ###
  hybrid.pvals <- c()
  for (i in seq_len(length(point_estimates))) {
    lambda <- sort(eigenvalues[, i], decreasing = T)
    if (all(lambda == 0)) {
      hybrid.pvals[i] <- NA
      next
    }
    Davies_Method <- suppressWarnings(davies(point_estimates[i], lambda = lambda, acc = acc, lim = 1e6))
    if (Davies_Method$ifault == 0) {
      hybrid.pvals[i] <- 2 * min(Davies_Method$Qq, 1 - Davies_Method$Qq)
    } else {
      warning(paste("Davies function exited with error code: ", Davies_Method$ifault))
      warning("Using saddlepoint approximation.")
      saddle <- saddlepoint_approximation(point_estimates[i], lambda = lambda)
      hybrid.pvals[i] <- 2 * min(saddle, 1 - saddle)
    }
  }
  return(hybrid.pvals)
}

# saddlepoint approx: adapted from https://github.com/cran/survey v4.0 and https://github.com/baolinwu/mkatr v0.1.0
saddlepoint_approximation = function(x, lambda) {
  d = max(lambda)
  lambda = lambda / d
  x = x / d
  k0 = function(zeta) - sum(log(1 - 2 * zeta * lambda)) / 2
  kprime0 = function(zeta) sapply(zeta, function(zz) sum(lambda / (1 - 2 * zz * lambda)))
  kpprime0 = function(zeta) 2*sum(lambda^2 / (1 - 2 * zeta * lambda)^2)
  n = length(lambda)
  if (any(lambda < 0)) {
    lmin = max(1 / (2 * lambda[lambda < 0])) * 0.99999
  } else if (x > sum(lambda)) {
    lmin = -0.01
  } else {
    lmin = -length(lambda) / (2 * x)
  }
  lmax = min(1 / (2 * lambda[lambda > 0])) * 0.99999
  hatzeta = uniroot(function(zeta) kprime0(zeta) - x, lower = lmin, upper = lmax, tol = 1e-08)$root
  w = sign(hatzeta) * sqrt(2 * (hatzeta * x - k0(hatzeta)))
  v = hatzeta * sqrt(kpprime0(hatzeta))
  if(abs(hatzeta) < 1e-4){
    return(NA)
  } else {
    return( pnorm(w + log(v / w) / w, lower.tail=FALSE) )
  }
}

sumlog <- function(pvalues) {
    df <- 2 * length(pvalues)
    fisherp <- pchisq(-2 * sum(log(pvalues)), df, lower.tail = FALSE)
    return(fisherp)
}
