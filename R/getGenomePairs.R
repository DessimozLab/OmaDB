#' Retrieves the pairwise relations among two genomes
#'
#' This function retrieves the pairwise relations among two genomes from
#' the OMA Browser database. The relations are orthologs in case the 
#' genomes are different and "close paralogs" and "homoeologs" in case 
#' they are the same.
#' 
#' By using the parameters chr1 and chr2, one can limit the relations 
#' to a certain chromosome for one or both genomes. The id of the 
#' chromosome corresponds to the chromosome ids from the 
#' \code{\link{getGenome}} result.
#' 
#' The rel_type parameter further limits the returned relations to a
#' specific subtype of orthologs (i.e. "1:1", "1:n", "m:1", "m:n") or 
#' - within a genome to either "close paralogs" or "homeologs".
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
#' getGenomePairs(genome_id1='YEAST',genome_id2='ASHGO')

getGenomePairs <- function(genome_id1, genome_id2, chr1 = NULL, chr2 = NULL, rel_type = NULL,
    ...) {
    if (missing(genome_id1) || missing(genome_id2)) {
        stop("You must provide IDs for both genomes.")
    }
    if (!is.null(rel_type)) {
        if (!(rel_type %in% c("1:1", "1:n", "m:n", "m:1", "close paralog", "homeolog"))) {
            stop("You must enter a valid relationship type.")
        }
    }

    url <- urlGenerator(endpoint = "pairs", id = genome_id1, detail = genome_id2,
        chr1 = chr1, chr2 = chr2, rel_type = rel_type)

    return(requestFactory(url, ...))
}


