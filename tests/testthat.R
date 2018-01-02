library(testthat)
library(roma)


test_that("Testing the correct format returned", {
    expect_equal(class(getData("group","YEAST58")$members), "data.frame")
    expect_equal(class(getXref(pattern="MAL")), "data.frame")
    expect_equal(class(getData("genome","YEAST")), "list")
    expect_equal(class(getGenomeAlignment("ASHGO","YEAST")), "data.frame")

	})
 	

	
test_that("Testing the URL validity", {
  expect_equal(status_code(httr::GET("https://omabrowser.org/api")), 200 ) #checking that the server is live and accesible

})



