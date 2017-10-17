#' Get Whole Genome Alignment Function
#' 
#' The function to obtain the list of orthologs for 2 whole genomes.
#'
#' @param genome_id1 an identifier for the first genome, which can be either its taxon id or UniProt species code
#' @param genome_id2 an an identifier for the second genome, which can be either its taxon id or UniProt species code
#' @param chr1 the chromosome of interest for the first genome 
#' @param chr2 the chromosome of interest for the second genome 
#' @return a dataframe containing information about both the entries in the orthologous pair and their relationship
#' @export
#' @examples
#' getGenomeAlignment(genome_id1="YEAST",genome_id2="ASHGO")
#' getGenomeAlignment(genome_id1="YEAST",genome_id2="ASHGO",chr1="1")


getGenomeAlignment <- function(genome_id1,genome_id2,chr1=NULL,chr2=NULL){
	if(missing(genome_id1) || missing(genome_id2)){
		stop("You must provide IDs for both genomes.")
	}
	if(!is.null(chr1) && !is.null(chr2)) {
		
		url = paste0(API_URL,"/pairs/",genome_id1,"/",genome_id2,"/?chr1=",chr1,"&chr2=",chr2)
	}
	if(!is.null(chr1) && is.null(chr2)) {
		
		url = paste0(API_URL,"/pairs/",genome_id1,"/",genome_id2,"/?chr1=",chr1)
	}
	if(is.null(chr1) && !is.null(chr2)) {
		
		url = paste0(API_URL,"/pairs/",genome_id1,"/",genome_id2,"/?chr2=",chr2)
	}
	else{
		url = paste0(API_URL,"/pairs/",genome_id1,"/",genome_id2,"/")
	}
	
	return(requestFactory(url))	
}