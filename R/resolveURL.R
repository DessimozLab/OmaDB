#' Get the Further Information behind the URL Function
#' 
#' The function to obtain further information from a given url. 
#'
#' @param url_field the url of interest
#' @param pagination boolean for which the default is TRUE
#' @return an object containing JSON keys as attributes 
#' @export
#' @examples
#' resolveURL(url_field="http://omadev.cs.ucl.ac.uk/api/protein/YEAST58/ontology/")

resolveURL <- function (url_field,pagination=TRUE) {
	return(requestFactory(url_field,pagination))
}