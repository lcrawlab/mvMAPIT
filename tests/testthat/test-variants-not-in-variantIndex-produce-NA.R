test_that(
    "pvalues and pve shows NA for variants not in variantIndex", {
        # given
        p <- 3
        n <- 10
        d <- 3
        set.seed(853)
        X <- matrix(
            runif(p * n),
            ncol = p
        )
        Y <- matrix(
            runif(d * n),
            ncol = d
        )
        variantIndex <- c(1, 3)
        otherIndex <- 2
        # when
        mapit <- MvMAPIT(
            t(X),
            t(Y),
            test = "normal", cores = 1, variantIndex = variantIndex, logLevel = "DEBUG"
        )
        # then
        otherp <- mapit$pvalues %>%
            filter(id == otherIndex)
        otherpve <- mapit$pves %>%
            filter(id == otherIndex)
        expect_true(all(is.na(otherp[, "p", drop = TRUE])))
        expect_true(all(is.na(otherpve[,"PVE", drop = TRUE])))
    }
)

test_that(
    "pve shows NA for covariance interactions", {
        # given
        p <- 3
        n <- 10
        d <- 3
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
        # then
        covint <- mapit$pves %>%
            filter(trait %in% c("P2*P1", "P3*P1", "P3*P2"))
        expect_true(all(is.na(covint[, "PVE", drop = TRUE])))
    }
)

