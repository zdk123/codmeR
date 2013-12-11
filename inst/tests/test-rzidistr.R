test_that("rzipois length equals input", {
    Ns <- seq(1, 1000,10)
    lapply(Ns, function(n) expect_equal(n, length(rzipois(n, 5, 0))))
})

test_that("rmvzipois dimensions match", {
    Ns <- seq(2, 1000, 200)
    Ds <- seq(2, 100, 10)
    for (n in Ns) {
        for (d in Ds) {
            data <- rmvzipois(n, mu=rep(20, d), diag(d))
            expect_equal(dim(data), c(n, d))
        }
    }

})
