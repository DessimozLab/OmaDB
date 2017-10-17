#' Get the Further Information behind the URL Function
#' 
#' The function to obtain further information from a given url. 
#'
#' @param url_field the url of interest
#' @return a data.frame containing the information behind an URL
#' @export
#' @examples
#' resolveURL(url_field="http://omadev.cs.ucl.ac.uk/api/protein/YEAST58/ontology/")

resolveURL <- function (url_field) {
	return(requestFactory(url_field))
}