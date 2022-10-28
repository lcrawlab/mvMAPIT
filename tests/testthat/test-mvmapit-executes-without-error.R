test_that(
    "test = 'normal'. ", {
        # given
        p <- 2
        n <- 10
        d <- 3
        pvalues <- tidyr::tibble(
           id = rep(as.character(c(1:p)), each = 6),
           trait = rep(c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3"), p),
           p =  c(0.4990573, 0.4478648, 0.9574136, 0.4662016, 0.4782672,
                  0.1381317, 0.5015375, 0.4619467, 0.1347061,
                  0.450717, 0.640529, 0.2410251)
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
        mapit <- mvmapit(
            t(X),
            t(Y),
            test = "normal", cores = 1, logLevel = "DEBUG"
        )
        # then
        expect_equal(mapit$pvalues, pvalues, tolerance = 1e-04)
    }
)

test_that(
    "test = davies. ", {
        # given
        p <- 2
        n <- 10
        d <- 3
        pvalues <- tidyr::tibble(
           id = rep(as.character(c(1:p)), each = 6),
           trait = rep(c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3"), p),
           p = c(0.01624319, NA, 0.6531582, NA, NA, 0.419213,
                 0.02694842, NA, 0.3977345, NA, NA, 0.490887)
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
        mapit <- mvmapit(
            t(X),
            t(Y),
            test = "davies", cores = 1, logLevel = "INFO"
        )
        # then
        expect_equal(mapit$pvalues, pvalues, tolerance = 1e-04)
    }
)

test_that(
    "test = hybrid", {
        # given
        p <- 2
        n <- 10
        d <- 3
        pvalues <- tidyr::tibble(
           id = rep(as.character(c(1:p)), each = 6),
           trait = rep(c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3"), p),
           p =  c(0.4990573, 0.4478648, 0.9574136, 0.4662016, 0.4782672,
                  0.1381317, 0.5015375, 0.4619467, 0.1347061,
                  0.450717, 0.640529, 0.2410251)
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
        mapit <- mvmapit(
            t(X),
            t(Y),
            test = "hybrid", cores = 1, logLevel = "INFO"
        )
        # then
        expect_equal(mapit$pvalues, pvalues, tolerance = 1e-04)
    }
)

test_that(
    "C is not NULL. ", {
        # given
        p <- 4
        n <- 10
        d <- 3
        pvalues <- tidyr::tibble(
           id = rep(as.character(c(1:p)), each = 6),
           trait = rep(c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3"), p),
           p =  c(
                0.6876487, 0.2148062, 0.5640931, 0.1657485, 0.2837563, 0.5020969,
                0.8920097, 0.9107812, 0.9787608, 0.6248188, 0.275113, 0.4994958,
                0.5868067, 0.5128342, 0.3823134, 0.874728, 0.2352273, 0.688964,
                0.3184337, 0.5047131, 0.6774045, 0.307193, 0.8257162, 0.5527816
           )
        )
        set.seed(29)
        X <- matrix(
            runif(p * n),
            ncol = p
        )
        Y <- matrix(
            runif(d * n),
            ncol = d
        )
        C <- matrix(
            runif(n * n),
            ncol = n
        )
        # when
        mapit <- mvmapit(
            t(X),
            t(Y),
            C = C, test = "hybrid", accuracy = 1e-05, cores = 1, logLevel = "ERROR"
        )
        # then
        expect_equal(mapit$pvalues, pvalues, tolerance = 1e-04)
    }
)

test_that(
    "test = 'normal', C is not NULL. ", {
        # given
        p <- 4
        n <- 10
        d <- 3
        pvalues <- tidyr::tibble(
           id = rep(as.character(c(1:p)), each = 6),
           trait = rep(c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3"), p),
           p =  c(
                0.6876487, 0.2148062, 0.5640931, 0.1657485, 0.2837563, 0.5020969,
                0.8920097, 0.9107812, 0.9787608, 0.6248188, 0.275113, 0.4994958,
                0.5868067, 0.5128342, 0.3823134, 0.874728, 0.2352273, 0.688964,
                0.3184337, 0.5047131, 0.6774045, 0.307193, 0.8257162, 0.5527816
           )
        )
        set.seed(29)
        X <- matrix(
            runif(p * n),
            ncol = p
        )
        Y <- matrix(
            runif(d * n),
            ncol = d
        )
        C <- matrix(
            runif(n * n),
            ncol = n
        )
        # when
        mapit <- mvmapit(
            t(X),
            t(Y),
            C = C, test = "normal", cores = 1, logLevel = "ERROR"
        )
        # then
        expect_equal(mapit$pvalues, pvalues, tolerance = 1e-04)
    }
)

test_that(
    "C is not NULL, test = 'davies'. ", {
        # given
        p <- 4
        n <- 10
        d <- 3
        pvalues <- tidyr::tibble(
           id = rep(as.character(c(1:p)), each = 6),
           trait = rep(c("P1*P1", "P2*P1", "P2*P2", "P3*P1", "P3*P2", "P3*P3"), p),
           p =  c(
                0.6977266, NA,0.4000627, NA, NA,0.3035032,
                0.4784944, NA,0.2936015, NA, NA,0.4665585,
                0.9274069, NA,0.1260672, NA, NA,0.1060084,
                0.4216691, NA,0.1615168, NA, NA,0.1526920
           )
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
        C <- matrix(
            runif(n * n),
            ncol = n
        )
        # when
        mapit <- mvmapit(
            t(X),
            t(Y),
            C = C, test = "davies", accuracy = 1e-05, cores = 1, logLevel = "ERROR"
        )
        # then
        expect_equal(mapit$pvalues, pvalues, tolerance = 1e-04)
    }
)

test_that(
    "test = 'davies'. , d = 1", {
        # given
        p <- 10
        n <- 4
        d <- 1
        pvalues <- tidyr::tibble(
           id = as.character(c(1:p)),
           trait = rep("P1", p),
           p =  c(
                0.7080633,
                0.0000000,
                0.0000000,
                0.1740376,
                0.2393295,
                0.2031174,
                0.0000000,
                0.1372735,
                0.1017353,
                0.2053481
           )
        )
        set.seed(20)
        X <- matrix(
            runif(p * n),
            ncol = p
        )
        Y <- matrix(
            runif(d * n),
            ncol = d
        )
        C <- matrix(
            runif(n * n),
            ncol = n
        )
        # when
        mapit <- mvmapit(
            t(X),
            t(Y),
            C = C, test = "davies", accuracy = 1e-05, cores = 1, logLevel = "ERROR"
        )
        # then
        expect_equal(mapit$pvalues, pvalues, tolerance = 1e-04)
    }
)

test_that(
    "test = 'hybrid'., d = 1 ", {
        # given
        p <- 2
        n <- 10
        d <- 1
        pvalues <- tidyr::tibble(
           id = rep(as.character(c(1:p))),
           trait = rep(c("P1"), p),
           p =  c(0.499, 0.502 )
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
        mapit <- mvmapit(
            t(X),
            t(Y),
            test = "hybrid", cores = 1, logLevel = "DEBUG"
        )
        # then
        print((mapit$pvalues))
        expect_equal(mapit$pvalues, pvalues, tolerance = 1e-03)
    }
)

