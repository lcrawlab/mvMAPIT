test_that(
    "test that all pvalues equal if phenotypes equal", {
        # given
        p <- 20
        n <- 10
        set.seed(853)
        X <- matrix(
            runif(p * n),
            ncol = p
        )
        y <- runif(n)
        Y <- as.matrix(cbind(y, y))
        # when
        mapit <- MvMAPIT(
            t(X),
            t(Y),
            test = "normal", cores = 1, logLevel = "DEBUG"
        )
        # then
        expect_equal(mapit$pvalues[, 1], mapit$pvalues[, 2], tolerance = 1e-08)
        expect_equal(mapit$pvalues[, 1], mapit$pvalues[, 3], tolerance = 1e-08)
    }
)

