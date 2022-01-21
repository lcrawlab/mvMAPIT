library('mvMAPIT')
# given
original <- readRDS("original_MAPIT.rds")
X <- original$genotype
Y <- original$phenotype
normal.pvalues <- original$normal$pvalues
davies.pvalues <- original$davies$pvalues
tolerance <- 0.0001
# when
Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); X=t((t(X)-Xmean)/Xsd)
mapit.normal <- MvMAPIT(t(X),
                 (Y),
                 test = 'normal',
                 cores = 4,
                 phenotypeCovariance = 'combinatorial',
                 logLevel = "INFO")
mapit.davies <- MvMAPIT(t(X),
                 t(Y),
                 test = 'davies',
                 cores = 4,
                 phenotypeCovariance = 'combinatorial',
                 logLevel = "INFO")
# then
normal <- as.vector(mapit.normal$pvalues)
names(normal.pvalues) <- NULL
names(normal) <- NULL
normal.all.equal <- all.equal(normal.pvalues, normal)
davies <- as.vector(mapit.davies$pvalues)
names(davies.pvalues) <- NULL
names(davies) <- NULL
davies.all.equal <- all.equal(davies.pvalues, davies)

normal.diff.counter <- 0
normal.na.counter <- c()
normal.diff.index <- c()

# comparison necessary because all.equal returns string if FALSE
if(normal.all.equal == TRUE) {
    message("Normal method p-values are equal for MAPIT and mvMAPIT.")
} else {
    warning("Normal method p-values differ between MAPIT and mvMAPIT.")
    NORMAL <- cbind(normal.pvalues, normal)
   for (i in seq_len(nrow(NORMAL))) {
        if(any(is.na(NORMAL[i, ]))) {
            normal.na.counter <- rbind(normal.na.counter, NORMAL[i, ])
            normal.diff.index <- c(normal.diff.index, i)
        } else if(abs(NORMAL[i, 1] - NORMAL[i, 2]) > tolerance) {
            print(NORMAL[i,])
            normal.diff.counter <- normal.diff.counter + 1
            normal.diff.index <- c(normal.diff.index, i)
        }
    }
}

davies.diff.counter <- 0
davies.na.counter <- c()
davies.diff.index <- c()

if(davies.all.equal == TRUE) {
    message("Davies method p-values are equal for MAPIT and mvMAPIT.")
} else {
    warning("Davies method p-values differ between MAPIT and mvMAPIT.")
    DAVIES <- cbind(davies.pvalues, davies)
    for (i in seq_len(nrow(DAVIES))) {
        if(any(is.na(DAVIES[i, ]))) {
            davies.na.counter <- rbind(davies.na.counter, DAVIES[i, ])
            davies.diff.index <- c(davies.diff.index, i)
        } else if (abs(DAVIES[i, 1] - DAVIES[i, 2]) > tolerance) {
            print(DAVIES[i,])
            davies.diff.counter <- davies.diff.counter + 1
            davies.diff.index <- c(davies.diff.index, i)
        }
    }
}

print(paste("Differences in normal:", normal.diff.counter))
print(paste("NA in normal:", nrow(normal.na.counter)))
print(paste("Differences in davies:", davies.diff.counter))
print(paste("NA in davies:", nrow(davies.na.counter)))

