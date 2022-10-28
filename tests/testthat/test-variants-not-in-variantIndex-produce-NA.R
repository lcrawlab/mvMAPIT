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
        mapit <- mvmapit(
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

