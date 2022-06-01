test_that(
    "Simulate multiple phenotypes returns apropriate phenotype object", {
        # given
        p <- 20
        f <- 10
        g <- 4
        n <- 5
        d <- 3
        X <- matrix(
            runif(p * n),
            ncol = p
        )

        # when
        data <- simulate_phenotypes(
            X, n_causal = f, n_trait_specific = g, n_pleiotropic = g, d = d, maf_threshold = 0,
            logLevel = "ERROR"
        )

        # then
        expect_equal(
            nrow(data$phenotype),
            n
        )
        expect_equal(
            ncol(data$phenotype),
            d
        )
        expect_equal(
            sum(is.na(data$phenotype)),
            0
        )  # no NA values in phenotype
        expect_equal(
            length(data),
            7
        )
    }
)

test_that(
    "Simulate multiple phenotypes returns causal SNPs", {
        # given
        p <- 20
        f <- 10
        g <- 4
        n <- 5
        d <- 3
        X <- matrix(
            runif(p * n),
            ncol = p
        )
        total_causal <- (f)  # single trait SNPs plus pleiotropic SNPs

        # when
        data <- simulate_phenotypes(
            X, n_causal = f, n_trait_specific = g, n_pleiotropic = g, d = d, maf_threshold = 0,
            logLevel = "ERROR"
        )

        # then
        expect_equal(
            nrow(data$additive),
            total_causal * d
        )
    }
)

test_that(
    "Simulate multiple phenotypes remove SNPs with low maf", {
        # given
        p <- 20
        f <- 10
        g <- 4
        n <- 5
        d <- 3
        maf <- 0.05 + 0.45 * runif(p)
        set.seed(12345)
        X <- matrix(
            (runif(p * n) >=
                maf) + (runif(p * n) >=
                maf), ncol = p
        )
        X[, 1] <- 0
        X[, 13] <- 0
        total_causal <- (f * p)  # single trait SNPs plus pleiotropic SNPs

        # when
        data <- simulate_phenotypes(
            X, n_causal = f, n_trait_specific = g, n_pleiotropic = g, d = d, maf_threshold = 0.05,
            logLevel = "ERROR"
        )
        data$snps.filtered

        # then
        expect_true(!(1 %in% data$snps.filtered) & !(13 %in% data$snps.filtered))
    }
)

test_that(
    "Simulation of groups with given ratio works as expected.", {
        # given
        p <- 50
        f <- 30
        g <- 8
        h <- 10
        n <- 5
        d <- 3
        X <- matrix(
            runif(p * n),
            ncol = p
        )
        total_causal <- (f)  # single trait SNPs plus pleiotropic SNPs
        correct_1 <- 2
        correct_2 <- 6
        correct_3 <- 2
        correct_4 <- 8

        # when
        data <- simulate_phenotypes(
            X, n_causal = f, n_trait_specific = g, n_pleiotropic = h, d = d, group_ratio_trait = 3,
            group_ratio_pleiotropic = 4, maf_threshold = 0, logLevel = "ERROR"
        )

        # then
        trait1 <- select(data$interactions, c("trait", "group1", "group2")) %>%
            filter(trait == 1) %>%
            ungroup() %>%
            tidyr::pivot_longer(col = c("group1", "group2")) %>%
            distinct()
        pleiotropic <- data$epistatic %>%
            filter(pleiotropic)
        t1g1 <- trait1 %>% filter(name == "group1")
        t1g2 <- trait1 %>% filter(name == "group2")
        result_1 <- length(setdiff(t1g1$value, pleiotropic$id))
        result_2 <- length(setdiff(t1g2$value, pleiotropic$id))
        result_3 <- length(intersect(t1g1$value, pleiotropic$id))
        result_4 <- length(intersect(t1g2$value, pleiotropic$id))
        expect_equal(result_1, correct_1)
        expect_equal(result_2, correct_2)
        expect_equal(result_3, correct_3)
        expect_equal(result_4, correct_4)
    }
)

test_that(
    "simulate_phenotypes can handle zero size n_trait_specific groups.", {
        # given
        p <- 50
        f <- 30
        g <- 0
        h <- 10
        n <- 5
        d <- 3
        X <- matrix(
            runif(p * n),
            ncol = p
        )
        total_causal <- (f)  # single trait SNPs plus pleiotropic SNPs
        correct_1 <- 0
        correct_2 <- 0

        # when
        data <- simulate_phenotypes(
            X, n_causal = f, n_trait_specific = g, n_pleiotropic = h, d = d, group_ratio_trait = 3,
            group_ratio_pleiotropic = 4, maf_threshold = 0, logLevel = "ERROR"
        )

        # then
        trait1 <- select(data$interactions, c("trait", "group1", "group2")) %>%
            filter(trait == 1) %>%
            ungroup() %>%
            tidyr::pivot_longer(col = c("group1", "group2")) %>%
            distinct()
        pleiotropic <- data$epistatic %>%
            filter(pleiotropic)
        t1g1 <- trait1 %>% filter(name == "group1")
        t1g2 <- trait1 %>% filter(name == "group2")
        result_1 <- length(setdiff(t1g1$value, pleiotropic$id))
        result_2 <- length(setdiff(t1g2$value, pleiotropic$id))
        expect_equal(result_1, correct_1)
        expect_equal(result_2, correct_2)
    }
)

test_that(
    "simulate_phenotypes can handle zero size n_pleiotropic groups.", {
        # given
        p <- 50
        f <- 30
        g <- 10
        h <- 0
        n <- 5
        d <- 3
        X <- matrix(
            runif(p * n),
            ncol = p
        )
        total_causal <- (f)  # single trait SNPs plus pleiotropic SNPs
        correct_1 <- 0
        correct_2 <- 0

        # when
        data <- simulate_phenotypes(
            X, n_causal = f, n_trait_specific = g, n_pleiotropic = h, d = d, group_ratio_trait = 3,
            group_ratio_pleiotropic = 4, maf_threshold = 0, logLevel = "ERROR"
        )
        pleiotropic <- data$epistatic %>%
            filter(pleiotropic)

        # then
        trait1 <- select(data$interactions, c("trait", "group1", "group2")) %>%
            filter(trait == 1) %>%
            ungroup() %>%
            tidyr::pivot_longer(col = c("group1", "group2")) %>%
            distinct()
        t1g1 <- trait1 %>% filter(name == "group1")
        t1g2 <- trait1 %>% filter(name == "group2")
        result_1 <- length(intersect(t1g1$value, pleiotropic$id))
        result_2 <- length(intersect(t1g2$value, pleiotropic$id))
        expect_equal(result_1, correct_1)
        expect_equal(result_2, correct_2)
    }
)

test_that(
    "simulate_phenotypes can handle zero size n_pleiotropic AND n_trait_specific groups.",
    {
        # given
        p <- 50
        f <- 30
        g <- 0
        h <- 0
        n <- 5
        d <- 3
        X <- matrix(
            runif(p * n),
            ncol = p
        )
        total_causal <- (f)  # single trait SNPs plus pleiotropic SNPs
        correct_1 <- NULL
        correct_2 <- NULL

        # when
        data <- simulate_phenotypes(
            X, n_causal = f, n_trait_specific = g, n_pleiotropic = h, d = d, group_ratio_trait = 3,
            group_ratio_pleiotropic = 4, maf_threshold = 0, logLevel = "ERROR"
        )

        # then
        expect_true(is.null(data$interactions))
        expect_true(is.null(data$epistatic))
    }
)

test_that(
    "test vector of heritability H2 <- c(0.6, 0.3)", {
        ind <- 100
        nsnp <- 100
        H2 <- c(0.6, 0.3)
        rho <- 0.5
        maf <- 0.05 + 0.45 * runif(nsnp)
        X <- (runif(ind * nsnp) <
            maf) + (runif(ind * nsnp) <
            maf)
        X <- matrix(
            as.double(X),
            ind, nsnp, byrow = TRUE
        )
        s <- 95345  # sample.int(10000, 1)
        sim <- simulate_phenotypes(
            X, n_causal = 30, n_pleiotropic = 6, n_trait_specific = 4, epistatic_correlation = 0.8,
            H2 = H2, rho = rho, logLevel = "ERROR", seed = s
        )
        expect_equal(length(sim$parameters$value), 16)

    }
)

test_that(
    "test run", {
        ind <- 100
        nsnp <- 100
        H2 <- 0.6
        rho <- 0.5
        maf <- 0.05 + 0.45 * runif(nsnp)
        X <- (runif(ind * nsnp) <
            maf) + (runif(ind * nsnp) <
            maf)
        X <- matrix(
            as.double(X),
            ind, nsnp, byrow = TRUE
        )
        s <- 95345  # sample.int(10000, 1)
        sim <- simulate_phenotypes(
            X, n_causal = 30, n_pleiotropic = 6, n_trait_specific = 4, epistatic_correlation = 0.8,
            H2 = H2, rho = rho, logLevel = "ERROR", seed = s
        )
    }
)
