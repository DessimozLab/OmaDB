library(testthat)

test_that("Testing the URL validity", {
  expect_equal(status_code(httr::GET("https://omabrowser.org/api")), 200 ) #checking that the server is live and accesible

})



