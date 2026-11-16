#' Retrieves the pairwise relations among two genomes
#'
#' This function retrieves the pairwise relations among two genomes from
#' the OMA Browser database. The relations are orthologs in case the
#' genomes are different and "close paralogs" and "homeologs" in case
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
#' @param genome_id1 an identifier for the first genome, which can be either
#'     its taxon id or UniProt species code
#' @param genome_id2 an identifier for the second genome, which can be
#'     either its taxon id or UniProt species code
#' @param chr1 the chromosome of interest for the first genome
#' @param chr2 the chromosome of interest for the second genome
#' @param rel_type the pairs relationship type
#' @param ... ignored; kept for backwards compatibility.
#' @return a data.frame containing information about both the entries in
#'     the pair (columns starting with \code{entry_1.} and \code{entry_2.})
#'     and their relationship.
#' @export
#' @examples
#' pairs <- getGenomePairs(genome_id1 = "YEAST", genome_id2 = "ASHGO",
#'     chr2 = "I", rel_type = "1:1")
#' head(pairs[, c("entry_1.omaid", "entry_2.omaid", "rel_type")])
getGenomePairs <- function(genome_id1, genome_id2, chr1 = NULL, chr2 = NULL,
        rel_type = NULL, ...) {
    if (missing(genome_id1) || missing(genome_id2)) {
        stop("You must provide IDs for both genomes.")
    }
    .checkScalarId(genome_id1)
    .checkScalarId(genome_id2)
    .checkRelType(rel_type,
        c("1:1", "1:n", "m:1", "m:n", "close paralog", "homeolog"))
    url <- .apiUrl("pairs", genome_id1, genome_id2,
        query = list(chr1 = chr1, chr2 = chr2, rel_type = rel_type))
    .toR(.omaGetList(url))
}
