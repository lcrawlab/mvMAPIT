test_that(
    "Fisher's combine method API", {
        # given
        p <- 2
        n <- 10
        d <- 3
        pvalues <- tidyr::tibble(
           id = rep(as.character(c(1:p)), each = 6),
           trait = rep(c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3"), p),
           p =  c(0.4990573, 0.4478648, 0.9574136, 0.4662016, 0.4782672,
                  0.1381317,  0.5015375, 0.4619467, 0.1347061,
                  0.450717, 0.640529, 0.2410251)
        )
        fisher <- tidyr::tibble(
            id = c("1", "2"),
            trait = c("fisher", "fisher"),
            p = c(0.612, 0.425)
        )
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
        mapit <- MvMAPIT(
            t(X),
            t(Y),
            test = "normal", cores = 1, logLevel = "DEBUG"
        )
        combined <- fishers_combined(mapit$pvalues)
        # then
        expect_equal(combined, fisher, tolerance = 1e-03)
    }
)

test_that(
    "Harmonic mean p combine method API", {
        # given
        p <- 2
        n <- 10
        d <- 3
        pvalues <- tidyr::tibble(
           id = rep(as.character(c(1:p)), each = 6),
           trait = rep(c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3"), p),
           p =  c(0.4990573, 0.4478648, 0.9574136, 0.4662016, 0.4782672,
                  0.1381317,  0.5015375, 0.4619467, 0.1347061,
                  0.450717, 0.640529, 0.2410251)
        )
        harmonic <- tidyr::tibble(
            id = c("1", "2"),
            trait = c("harmonic", "harmonic"),
            p = c(0.358, 0.308)
        )
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
        mapit <- MvMAPIT(
            t(X),
            t(Y),
            test = "normal", cores = 1, logLevel = "DEBUG"
        )
        combined <- harmonic_combined(mapit$pvalues)
        # then
        expect_equal(combined, harmonic, tolerance = 1e-03)
    }
)
