#' Get the taxonomy of the genomes in the OMA Browser
#'
#' The function obtains the taxonomic tree of the genomes in the OMA
#' Browser, either as a newick string (that can for example be converted into
#' a \code{phylo} object with \code{\link{getTree}}) or as a nested
#' object.
#'
#' The taxonomy can be restricted to a subtree by specifying its root, or
#' to the taxonomy induced by a set of member genomes or taxonomic levels.
#' Only one of \code{root} and \code{members} can be specified.
#'
#' @param root optional parameter, the root of the subtree of interest, i.e.
#'     a taxonomic level, a taxon id or a UniProt species code.
#' @param members optional parameter, vector (or comma-separated string) of
#'     NCBI taxon ids, species codes or taxonomic levels that should be
#'     included in the induced taxonomy.
#' @param newick logical, whether to return the taxonomy in newick format
#'     (default) or as a nested object.
#' @param collapse logical, whether levels with a single child should be
#'     collapsed (default TRUE).
#' @return an \code{\link{omadb_obj}}. If \code{newick = TRUE}, it contains
#'     the attributes \code{newick} (labels containing special characters
#'     are quoted) and \code{root_taxon}. Otherwise it contains the attributes
#'     \code{name}, \code{id} and \code{children}, a data.frame with one row
#'     per child taxon, where a list column \code{children} again contains
#'     the children of internal nodes and the column \code{code} the species
#'     code of leaves.
#' @export
#' @examples
#' getTaxonomy(members = c("YEAST", "ASHGO"))$newick
#' tax <- getTaxonomy(root = "Alveolata", newick = FALSE)
#' tax$name
getTaxonomy <- function(root = NULL, members = NULL, newick = TRUE,
        collapse = TRUE) {
    if (!is.logical(newick) || length(newick) != 1L || is.na(newick)) {
        stop("newick parameter must be TRUE or FALSE.")
    }
    if (!is.null(root) && !is.null(members)) {
        stop("Only one of root and members can be specified.")
    }
    type <- if (newick) "newick" else "dictionary"
    query <- list(type = type, collapse = if (!collapse) FALSE,
        newick_quote_labels = if (newick) TRUE)
    if (!is.null(members)) {
        query$members <- paste(members, collapse = ",")
    }
    .toR(.omaGet(.apiUrl("taxonomy", root, query = query)))
}
