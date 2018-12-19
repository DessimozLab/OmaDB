#' Get Genome Pairs Function
#' 
#' The function to obtain the list of orthologs for 2 whole genomes.
#'
#' @param genome_id1 an identifier for the first genome, which can be either its taxon id or UniProt species code
#' @param genome_id2 an an identifier for the second genome, which can be either its taxon id or UniProt species code
#' @param chr1 the chromosome of interest for the first genome 
#' @param chr2 the chromosome of interest for the second genome 
#' @param rel_type the pairs relationship type
#' @param ... qwargs
#' @return a dataframe containing information about both the entries in the orthologous pair and their relationship
#' @export
#' @examples
#' getGenomePairs(genome_id1="YEAST",genome_id2="ASHGO")



getGenomePairs <- function(genome_id1,genome_id2,chr1=NULL,chr2=NULL,rel_type=NULL, ...){
	if(missing(genome_id1) || missing(genome_id2)){
		stop("You must provide IDs for both genomes.")
	}
	if(!is.null(rel_type)){
		if(!(rel_type %in% c("1:1","1:n","m:n","m:1"))){
			stop("You must enter a valid relationship type.")
			}	
	}
	
	url = urlGenerator(endpoint='pairs', id=genome_id1, detail=genome_id2, chr1=chr1, chr2=chr2, rel_type=rel_type)

	return(requestFactory(url, ... ))	
}


