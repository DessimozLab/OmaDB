#' Get the Data Function
#' 
#' The function to obtain the information available for a single entry in the datase.
#'
#' @param type the type for the entry to be returned - either protein, genome, group or hog
#' @param id an identifier for the entry to be returned. For more information, see the "Get started with Roma" vignette.
#' @return an object containing the JSON keys as attributes
#' @export
#' @examples
#' getData(type = "Protein", id="YEAST00001")
#' getData(type = "Group", id="YEAST00001")


getData <- function(type, id=NULL){
	if(missing(type)){
		stop("You must provide a valid object type.")
	}
	url = urlGenerator(type=type,id=id)
	
	check_response(url)

	return(requestFactory(url))
}