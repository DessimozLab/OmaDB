#' Get the Data Function
#' 
#' The function to obtain the information available for a single entry in the datase.
#'
#' @param type the type for the entry to be returned - either protein, genome, group or hog
#' @param id an identifier for the entry to be returned. For more information, see the "Get started with Roma" vignette.
#' @return an object containing the JSON keys as attributes
#' @export
#' @examples
#' getData(type = "protein", id="YEAST00001")
#' getData(type = "group", id="YEAST00001")


getData <- function(type, id=NULL){

	type = tolower(type)

	if(missing(type) || !(type %in% list("group","protein","genome"))){
		stop("You must provide a valid object type.")
	}

	if(missing(id) || grepl(",",id) || length(id)!= 1){
		stop("You must provide a valid object id.")
	}

		
	url = urlGenerator(type=type,id=id)
		
	return(requestFactory(url))
	

	
}