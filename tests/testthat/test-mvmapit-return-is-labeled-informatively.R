test_that(
    "pairwise test with d = 3. test = 'normal'.", {
        # given
        p <- 10
        n <- 5
        d <- 3
        variants <- sprintf("SNP%s", 1:p)
        phenotypes <- sprintf("Q%s", 1:d)
        traits <- c("Q1*Q1", "Q2*Q1", "Q2*Q2", "Q3*Q1", "Q3*Q2", "Q3*Q3")
        ids <- rep(as.character(variants), each = 6)
        trait_values <- rep(traits, p)
        pvalues <- tidyr::tibble(
           id = ids,
           trait = trait_values,
           p = rep(1, 60)
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
        colnames(X) <- variants
        colnames(Y) <- phenotypes
        # when
        mapit <- mvmapit(
            t(X),
            t(Y),
            test = "normal", cores = 1, logLevel = "ERROR"
        )
        # then
        expect_equal(
            mapit$pvalues[ ,"id", drop = TRUE],
            ids
        )
        expect_equal(
            mapit$pvalues[ ,"trait", drop = TRUE],
            trait_values
        )
    }
)

