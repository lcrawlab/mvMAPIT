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
        p1p1 <- mapit$pvalues %>%
            filter(trait == "y*y") %>%
            select(-trait)
        p2p1 <- mapit$pvalues %>%
            filter(trait == "y.1*y") %>%
            select(-trait)
        p2p2 <- mapit$pvalues %>%
            filter(trait == "y.1*y.1") %>%
            select(-trait)
        expect_equal(p1p1, p2p1, tolerance = 1e-08)
        expect_equal(p1p1, p2p2, tolerance = 1e-08)
    }
)

