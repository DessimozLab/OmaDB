
#' Get the Protein data function
#' 
#' The function to obtain the information available for a single protein or multiple proteins in the database.
#'
#' @param id identifier(s) for the entry or entries to be returned. a character string if single entry or a vector if multiple.
#' @param attribute an extra attribute to be returned.
#' @param matchPartially boolean. If set to TRUE, in case there is no exact match is found in the database it will return top hits.
#' @return an object containing the JSON keys as attributes or a dataframe
#' @export
#' @examples
#' getProtein(id="YEAST00001")
#' getProtein(id="YEAST00001",attribute='ontology')
#' getProtein(id=c("YEAST00001","YEAST00002","YEAST00012"))
#' getProtein(id=c("YEAST00001","YEAST00002","YEAST00012"),attribute='ontology')
#' getProtein(id="MAL", matchPartially=TRUE)


getProtein <- function(id, attribute = NULL, matchPartially = FALSE){

	if(!is.null(attribute) && !(attribute %in% c('domains','homeologs','ontology','orthologs'))){
		stop("You must provide a valid attribute.")
	}

	if(matchPartially==TRUE){
	    if (!length(id)==1){ stop("You can only search with a single partial ID"); }
		url = urlGenerator(endpoint="xref", search=id)
	}

	else if(length(id)==1){
		url = urlGenerator(endpoint='protein', id=id, detail=attribute)	
	}

	else if(length(id) > 1){
		body = jsonlite::toJSON(list(ids=id, auto_unbox=TRUE))
		url = urlGenerator(endpoint = "protein", id = "bulk_retrieve")
		data = requestFactory(url = url, body = body)
		names(data) <- id

		if(is.null(attribute)){
			return(data)
		} else {
			attribute_data = lapply(data, function(x) if(grepl('https://',x[[attribute]])){requestFactory(x[[attribute]]) } else{ x[[attribute]] })
			return(attribute_data)
		}
	}
	return(requestFactory(url))
}

#' Get the Genome data function
#' 
#' The function to obtain the information available for a single protein or multiple proteins in the datase.
#'
#' @param id identifier(s) for the entry or entries to be returned. a character string if single entry or a vector if multiple.
#' @param attribute an extra attribute to be returned (proteins)
#' @return an object containing the JSON keys as attributes or a dataframe 
#' @export
#' @examples
#' getGenome(id="HUMAN")
#' getGenome(id="HUMAN",attribute='proteins')

getGenome <- function(id,attribute=NULL){

	if(!is.null(attribute) && !(attribute %in% c('proteins'))){
		stop("You must provide a valid attribute.")
	}

	url = urlGenerator(endpoint='genome', id=id, detail=attribute)	

	return(requestFactory(url))
}

#' Get the OMA group data function
#' 
#' The function to obtain the information available for a single protein or multiple proteins in the datase.
#'
#' @param id identifier(s) for the entry or entries to be returned. a character string if single entry or a vector if multiple.
#' @param attribute an extra attribute to be returned (close_groups)
#' @return an object containing the JSON keys as attributes or a dataframe
#' @export
#' @examples
#' getOMAGroup(id="58")
#' getOMAGroup(id="58",attribute='close_groups')

getOMAGroup <- function(id, attribute=NULL){

	if(!is.null(attribute) && !(attribute %in% c('close_groups'))){
		stop("You must provide a valid attribute.")
	}
	url = urlGenerator(endpoint='group', id=id, detail=attribute)	
	return(requestFactory(url))
}

