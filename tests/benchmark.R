library(mvMAPIT)

N <- c(16000)
p <- 10
d <- 1
for(n in N) {
    set.seed(853)
    X <- matrix(
        runif(p * n),
        ncol = p
    )
    Y <- matrix(
        runif(d * n),
        ncol = d
    )
    # when
    mapit <- mvmapit(
        t(X),
        t(Y),
        test = "normal", cores = 1, logLevel = "DEBUG"
    )
    # then
    print(mapit$duration)
}
