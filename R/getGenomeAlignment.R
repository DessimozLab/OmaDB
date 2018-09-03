#' Get Whole Genome Alignment Function
#' 
#' The function to obtain the list of orthologs for 2 whole genomes.
#'
#' @param genome_id1 an identifier for the first genome, which can be either its taxon id or UniProt species code
#' @param genome_id2 an an identifier for the second genome, which can be either its taxon id or UniProt species code
#' @param chr1 the chromosome of interest for the first genome 
#' @param chr2 the chromosome of interest for the second genome 
#' @param rel_type the pairs relationship type
#' @param per_page the number of instances to be returned or 'all'. default is set to a 100. 
#' @return a dataframe containing information about both the entries in the orthologous pair and their relationship
#' @export
#' @examples
#' getGenomeAlignment(genome_id1="YEAST",genome_id2="ASHGO")
#' getGenomeAlignment(genome_id1="YEAST",genome_id2="ASHGO",chr1="1")


getGenomeAlignment <- function(genome_id1,genome_id2,chr1=NULL,chr2=NULL,per_page=NULL,rel_type=NULL){
	if(missing(genome_id1) || missing(genome_id2)){
		stop("You must provide IDs for both genomes.")
	}
	if(!is.null(rel_type)){
		if(!(rel_type %in% c("1:1","1:n"))){
			stop("You must enter a valid relationship type.")
			}	
	}

	if(!is.null(chr1) && !is.null(chr2) && !is.null(rel_type)) {
		
		url = urlGenerator(type = 'pairs', id=genome_id1, detail = genome_id2,query_param1='chr1',query_param1_value=chr1,query_param2='chr2',query_param2_value=chr2,query_param3='rel_type',query_param3_value=rel_type)
	}


	if(is.null(chr1) && is.null(chr2) && !is.null(rel_type)) {
		
		url = urlGenerator(type = 'pairs', id=genome_id1, detail = genome_id2,query_param1='rel_type',query_param1_value=rel_type)
	}

	if(!is.null(chr1) && !is.null(chr2) && is.null(rel_type)) {
		
		url = urlGenerator(type = 'pairs', id=genome_id1, detail = genome_id2,query_param1='chr1',query_param1_value=chr1,query_param2='chr2',query_param2_value=chr2)
	}


	if(!is.null(chr1) && is.null(chr2) && is.null(rel_type)) {
		
		url = urlGenerator(type = 'pairs', id=genome_id1, detail = genome_id2,query_param1='chr1',query_param1_value=chr1)
	}

	if(is.null(chr1) && !is.null(chr2) && is.null(rel_type)) {
		
		url = urlGenerator(type = 'pairs', id=genome_id1, detail = genome_id2,query_param1='chr2',query_param1_value=chr2)
	}

	if(is.null(chr1) && !is.null(chr2) && !is.null(rel_type)) {
		
		url = urlGenerator(type = 'pairs', id=genome_id1, detail = genome_id2,query_param1='chr2',query_param1_value=chr2,query_param3='rel_type',query_param3_value=rel_type)
	}

	if(!is.null(chr1) && is.null(chr2) && !is.null(rel_type)) {
		
		url = urlGenerator(type = 'pairs', id=genome_id1, detail = genome_id2,query_param1='chr1',query_param1_value=chr1,query_param3='rel_type',query_param3_value=rel_type)
	}

	if(is.null(chr1) && is.null(chr2) && is.null(rel_type)) {
		
		url = urlGenerator(type = 'pairs', id=genome_id1, detail = genome_id2)
	}
	
	return(requestFactory(url,per_page = per_page))	
}


