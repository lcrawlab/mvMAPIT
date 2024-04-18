test_that(
    "Toggle Projections Output Correlated", {
        # given
        p <- 2
        n <- 10
        d <- 1

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
            test = "normal", cores = 1, logLevel = "DEBUG", skipProjection = FALSE
        )
        
        mapit_no_projections <- mvmapit(
            t(X),
            t(Y),
            test = "normal", cores = 1, logLevel = "DEBUG", skipProjection = TRUE
        )
        print(mapit$pvalues$p)
        print( mapit_no_projections$pvalues$p)
        correlation <- cor(mapit$pvalues$p, mapit_no_projections$pvalues$p)
        print(correlation)
        
        # then expect_equal(correlation, smth, tolerance = 1e-03)
        
    }
)

#toggle projections off once and then not -- projections true or false
