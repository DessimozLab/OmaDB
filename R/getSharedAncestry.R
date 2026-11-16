#' Get the shared ancestry of two genomes
#'
#' The function retrieves summary statistics on the shared ancestry of two
#' genomes, i.e. the fraction of genes of the two genomes that have an
#' ortholog in the other genome, and the most recent common ancestor of the
#' two genomes.
#'
#' @param genome_id1 an identifier for the first genome, which can be
#'     either its taxon id or UniProt species code
#' @param genome_id2 an identifier for the second genome, which can be
#'     either its taxon id or UniProt species code
#' @return an \code{\link{omadb_obj}} with the attributes \code{fraction},
#'     \code{details} (a data.frame with the number of genes and orthologs
#'     per genome) and \code{mrca} (the most recent common ancestor).
#' @export
#' @examples
#' res <- getSharedAncestry("YEAST", "ASHGO")
#' res$fraction
#' res$details
getSharedAncestry <- function(genome_id1, genome_id2) {
    if (missing(genome_id1) || missing(genome_id2)) {
        stop("You must provide IDs for both genomes.")
    }
    .checkScalarId(genome_id1)
    .checkScalarId(genome_id2)
    url <- .apiUrl("summary", "shared_ancestry", genome_id1, genome_id2)
    .toR(.omaGet(url))
}
