#' Get loci for a given list of proteins
#'
#' Function to obtain the genomic loci of a set of proteins as a GRanges
#' object.
#'
#' @param proteins the proteins of interest. Either a character vector of
#'     protein ids, a data.frame with protein data (e.g. the members of an
#'     OMA group or HOG, which have the columns \code{omaid},
#'     \code{chromosome}, \code{locus.start}, \code{locus.end} and
#'     \code{locus.strand}), a list of such data.frames or an
#'     \code{\link{omadb_obj}} of a protein.
#' @return a \code{\link[GenomicRanges]{GRanges}} object with the loci of the
#'     proteins. The seqnames are the chromosomes, the names are the OMA ids
#'     of the proteins.
#' @export
#' @examples
#' loci <- getLocus(proteins = getOMAGroup("YEAST58")$members)
#' loci
#' getLocus(c("YEAST00001", "YEAST00002"))
getLocus <- function(proteins) {
    if (is.character(proteins)) {
        return(getProtein(id = proteins, attribute = "locus"))
    }
    if (is(proteins, "omadb_obj")) {
        return(.subset2(proteins, "locus"))
    }
    if (is.list(proteins) && !is.data.frame(proteins)) {
        if (!all(vapply(proteins, is.data.frame, logical(1)))) {
            stop("proteins must be ids, a data.frame or a list of ",
                "data.frames.")
        }
        proteins <- .rbindFill(proteins)
    }
    required <- c("omaid", "locus.start", "locus.end", "locus.strand")
    if (!is.data.frame(proteins) || !all(required %in% names(proteins))) {
        stop("proteins must contain the columns ",
            paste(required, collapse = ", "), ".")
    }
    chromosome <- proteins[["chromosome"]]
    if (is.null(chromosome)) {
        chromosome <- rep("unknown", nrow(proteins))
    }
    chromosome[is.na(chromosome)] <- "unknown"
    GenomicRanges::GRanges(
        seqnames = chromosome,
        ranges = IRanges::IRanges(start = as.integer(proteins$locus.start),
            end = as.integer(proteins$locus.end), names = proteins$omaid),
        strand = .toStrand(proteins$locus.strand))
}
