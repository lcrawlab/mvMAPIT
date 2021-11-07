library('mvMAPIT')
# given
original <- readRDS("original_MAPIT.rds")
X <- original$genotype
Y <- original$phenotype
normal.pvalues <- original$normal$pvalues
davies.pvalues <- original$davies$pvalues
# when
mapit.normal <- MvMAPIT(t(X),
                 (Y),
                 test = 'normal',
                 cores = 1,
                 phenotypeCovariance = 'combinatorial',
                 logLevel = "ERROR")
mapit.davies <- MvMAPIT(t(X),
                 t(Y),
                 test = 'davies',
                 cores = 1,
                 phenotypeCovariance = 'combinatorial',
                 logLevel = "ERROR")
# then
normal <- as.vector(mapit.normal$pvalues)
names(normal.pvalues) <- NULL
names(normal) <- NULL
all.equal(normal.pvalues, normal)
davies <- as.vector(mapit.davies$pvalues)
names(davies.pvalues) <- NULL
names(davies) <- NULL
all.equal(davies.pvalues, davies)
DAVIES <- cbind(davies.pvalues, davies)
NORMAL <- cbind(normal.pvalues, normal)
tolerance <- 0.0001
normal.diff.counter <- 0
normal.na.counter <- c()
for (i in seq_len(nrow(NORMAL))) {
    if(is.na(NORMAL[i, 1]) || is.na(NORMAL[i, 2])) {
        print(NORMAL[i, ])
        normal.na.counter <- rbind(normal.na.counter, NORMAL[i, ])
        next
    }
    if(abs(NORMAL[i, 1] - NORMAL[i, 2]) > tolerance)
        {
            print(NORMAL[i,])
            normal.diff.counter <- normal.diff.counter + 1
        }
}
davies.diff.counter <- 0
davies.na.counter <- c()
for (i in seq_len(nrow(DAVIES))) {
    if(is.na(DAVIES[i, 1]) || is.na(DAVIES[i, 2])) {
        print(DAVIES[i, ])
        davies.na.counter <- rbind(normal.na.counter, DAVIES[i, ])
        next
    }
    if(abs(DAVIES[i, 1] - DAVIES[i, 2]) > tolerance)
        {
            print(DAVIES[i,])
            davies.diff.counter <- normal.diff.counter + 1
        }
}
print(paste("Differences in normal:", normal.diff.counter))
print(paste("NA in normal:", nrow(normal.na.counter)))
print(paste("Differences in davies:", davies.diff.counter))
print(paste("NA in davies:", nrow(davies.na.counter)))

