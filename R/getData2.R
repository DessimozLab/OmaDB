#' Retrieve a protein from the OMA Browser
#' 
#' This function enables to retrieve information on one or several proteins from the 
#' OMA Browser database. 
#'
#' In its simplest form the function returns the base data of the query protein.
#' The query protein can be selected with any unique id, for example with a 
#' UniProtKB accession (P12345), an OMA id (YEAST00012), or a RefSeq id (NP_001226).
#' To retrieve more than one protein, you should pass a vector of IDs.
#'
#' Non-scalar properties of proteins such as their domains, GO annotations, 
#' orthologs or homeologs will get loaded upon accessing them, or if you only
#' need this information you can set the attribute parameter to the property name
#' and retrieve this information directly.
#'
#' @seealso For non-unique non-unique IDs or partial ID lookup, use [searchProtein()] instead.
#'
#' @param id Identifier(s) for the entry or entries to be returned. a character string if single entry or a vector if multiple.
#' @param attribute Instead of the protein, return the attribute property of the protein. Attriute needs to be one of 'domains', 'orthologs', 'ontology', or 'homoeologs'.
#' @return An object containing the JSON keys as attributes or a dataframe containing the non-scalar protein property.
#' @export
#' @examples
#' getProtein(id="YEAST00001")
#' getProtein(id="YEAST00001", attribute='orthologs')
#' getProtein(id=c("YEAST00001","YEAST00002","YEAST00012"))
#' getProtein(id=c("YEAST00001","YEAST00002","YEAST00012"), attribute='ontology')


getProtein <- function(id, attribute = NULL){

	if(!is.null(attribute) && !(attribute %in% c('domains','homeologs','ontology','orthologs'))){
		stop("You must provide a valid attribute.")
	}

	if(length(id)==1){
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

#' Retrieve an OMA Group from the OMA Browser
#' 
#' This function obtains an OMA Group from the OMA Browser database. An OMA Group
#' is defined to be a clique of proteins that are all orthologous to each other, 
#' i.e. they are all related through speciation events only. An OMA Group can thus
#' by definition not contain any inparalogs. It is a very stringent orthology
#' grouping approach.
#' OMA Groups are mostly useful to infer phylogenetic species tree where they
#' can be used as marker genes.
#'
#' Retrieving an OMA Group can be done using a group nr as id, its fingerprint
#' (a 7mer AA sequence which is unique to proteins in that group), a member 
#' protein id or any sequence pattern that is unique to the group.
#'
#' @param id An identifier for the group. See above for possible types of IDs.
#' @param attribute an extra attribute to be returned (close_groups)
#' @return an object containing the JSON keys as attributes or a dataframe
#' @export
#' @examples
#' getOMAGroup(id="58")
#' getOMAGroup(id="P12345")
#' getOMAGroup(id="NNRRGRI")
#' getOMAGroup(id="58", attribute='close_groups')

getOMAGroup <- function(id, attribute=NULL){

	if(!is.null(attribute) && !(attribute %in% c('close_groups'))){
		stop("You must provide a valid attribute.")
	}
	url = urlGenerator(endpoint='group', id=id, detail=attribute)	
	return(requestFactory(url))
}

