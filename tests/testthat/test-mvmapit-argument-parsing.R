test_that(
    "MvMapit can take a vector as phenotype input. test = normal", {
        # given
        p <- 10
        n <- 4
        pvalues <- tidyr::tibble(
           id = as.character(c(1:p)),
           trait = rep("P1", p),
           p =  rep(0.48001, p)
        )
        set.seed(853)
        X <- matrix(
            runif(p * n),
            ncol = p
        )
        Y <- c(runif(n))
        # when
        mapit <- mvmapit(
            t(X),
            Y, accuracy = 1e-05, cores = 1, logLevel = "ERROR"
        )
        # then
        expect_equal(mapit$pvalues, pvalues, tolerance = 0.01)
    }
)

test_that(
    "MvMapit can take a vector as phenotype input. test = davies", {
        # given
        p <- 10
        n <- 4
        pvalues <- tidyr::tibble(
           id = as.character(c(1:p)),
           trait = rep("P1", p),
           p =  rep(0.209, p)
        )
        set.seed(853)
        X <- matrix(
            runif(p * n),
            ncol = p
        )
        Y <- c(runif(n))
        # when
        mapit <- mvmapit(
            t(X),
            Y, test = "davies", accuracy = 1e-05, cores = 1, logLevel = "ERROR"
        )
        # then
        expect_equal(mapit$pvalues, pvalues, tolerance = 0.001)
    }
)
