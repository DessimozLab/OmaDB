#' Get the Data Function
#' 
#' The function to obtain the information available for a single entry in the datase.
#'
#' @param type the type for the entry to be returned - either protein, genome or group
#' @param id an identifier for the entry to be returned. For more information, see the "Get started with OmaDB" vignette.
#' @param attribute an extra attribute 
#' @return an object containing the JSON keys as attributes
#' @export
#' @examples
#' getData(type = "protein", id="YEAST00001")
#' getData(type = "group", id="YEAST00001")


getData <- function(type, id=NULL, attribute = NULL){

	type = tolower(type)

	if(missing(type) || !(type %in% list("group","protein","genome"))){
		stop("You must provide a valid object type.")
	}

	if(missing(id) || grepl(",",id) || length(id)!= 1){
		stop("You must provide a valid object id.")
	}

	if(is.null(attribute)){
		url = urlGenerator(type=type,id=id)

	}

	if(!is.null(attribute) && type!='protein'){
		stop('An attribute functionality is currently only available with type protein.')
	}

	if(!is.null(attribute) && !(attribute %in% c('domains','homeologs','ontology','orthologs'))){
		stop('You must choose a valid attribute.')
	}

	if(!is.null(attribute) && attribute %in% c('domains','homeologs','ontology')){
		url = urlGenerator(type=type,id = id, detail = attribute)
	}

	return(requestFactory(url))
	

	
}