test_that("running mean with constant x or position", {

  n <- 100
  x <- rnorm(n)
  pos <- rep(0, n)
  expect_equal( runningmean(pos, x, window=1), rep(mean(x), n) )

  mu <- mean(x)
  x <- rep(mu, n)
  pos <- runif(n, 0, 5)
  expect_equal( runningmean(pos, x, window=1), x)

})