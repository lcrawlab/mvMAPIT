test_that(
    "saddlepoint approximation and davies similar results", {
        # given
        x <- c(1, 5, 11, 15, 20)
        lambda <- rep(1, 10)
        # when
        davies_results <- c()
        saddle_results <- c()
        for (i in x) {
            saddle_results <- c(saddle_results, saddlepoint_approximation(i, lambda))
            davies_results <- c(davies_results, davies(i, lambda = lambda)$Qq)
        }
        # then
        expect_equal(saddle_results, davies_results, tolerance = 0.001)
    }
)
