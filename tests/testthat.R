library(testthat)
library(roma)
library(topGO)


test_that("Testing the correct format returned", {
    expect_equal(class(getData("group","YEAST58")$members), "data.frame")
    expect_equal(class(getXref(pattern="MAL")), "data.frame")
    expect_equal(class(getData("genome","YEAST")), "list")
    expect_equal(class(getGenomeAlignment("ASHGO","YEAST")), "data.frame")

   
    geneList = list(getData(type="protein",id="YEAST560"),getData(type="protein",id="YEAST346"))
    annotations = formatTopGO(geneList,format="geneID2GO")   
    expect_match(class(getTopGO(annotations=annotations, 
   								myInterestingGenes = list("YEAST00346"), 
   								format = "geneID2GO")), "topGOdata")

	})
 	

	
test_that("Testing the URL validity", {
  expect_equal(status_code(httr::GET("https://omabrowser.org/api")), 200 ) #checking that the server is live and accesible

})



