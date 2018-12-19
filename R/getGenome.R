#' Retrieve a genome from the OMA Browser database
#' 
#' This function obtains the basic information for one specific genome available
#' on the OMA Browser, or - if no id is provided - a dataframe with all available
#' genomes.
#' 
#' Ids can be either the scientific name of a species, the NCBI taxonomy id or 
#' the UniProtKB mnemonic species code.
#'
#' The optional argument attribute can be used to directly load the proteins 
#' belonging to the genome. Alternatively, you can access the proteins attribute
#' of the result which will transparently load the proteins from the OMA Browser.
#'
#' @param id A genome identifier. By default, all available genomes will be returned.
#' @param attribute An extra attribute to be returned (proteins)
#' @return an object containing the JSON keys as attributes or a dataframe 
#' @export
#' @examples
#' getGenome()
#' getGenome(id="HUMAN")
#' getGenome(id="HUMAN",attribute='proteins')

getGenome <- function(id=NULL, attribute=NULL){

	if(!is.null(attribute) && !(attribute %in% c('proteins'))){
		stop("You must provide a valid attribute.")
	}

	url = urlGenerator(endpoint='genome', id=id, detail=attribute)	

	return(requestFactory(url))
}