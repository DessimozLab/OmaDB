#' Get the Tree Object
#'
#' A convenience function to obtain a tree object from a newick tree,
#' essentially wraps \code{\link[ape]{read.tree}} from the ape package.
#'
#' @param newick The newick tree to be instantiated, or the result of
#'     \code{\link{getTaxonomy}} with \code{newick = TRUE}.
#' @return a tree object of class \code{phylo}
#' @export
#' @importFrom ape read.tree
#' @examples
#' taxonomy <- getTaxonomy(root = "Alveolata")
#' tree <- getTree(newick = taxonomy)
#' head(tree$tip.label)
getTree <- function(newick) {
    if (is.list(newick)) {
        newick <- .subset2(newick, "newick")
    }
    if (!is.character(newick) || length(newick) != 1L) {
        stop("newick must be a newick string or a taxonomy object.")
    }
    .readNewick(newick)
}

## Read a newick string with ape. ape does not support labels in double
## quotes (as returned by the OMA Browser for labels with special
## characters), so these are replaced by placeholders while parsing.
.readNewick <- function(newick) {
    m <- gregexpr('"[^"]*"', newick)
    labels <- regmatches(newick, m)[[1]]
    if (length(labels) == 0L) {
        return(ape::read.tree(text = newick))
    }
    keys <- sprintf("omadblabel%06d", seq_along(labels))
    regmatches(newick, m) <- list(keys)
    tree <- ape::read.tree(text = newick)
    lookup <- stats::setNames(substr(labels, 2L, nchar(labels) - 1L), keys)
    tree$tip.label <- .restoreLabels(tree$tip.label, lookup)
    if (!is.null(tree$node.label)) {
        tree$node.label <- .restoreLabels(tree$node.label, lookup)
    }
    tree
}

.restoreLabels <- function(x, lookup) {
    hit <- x %in% names(lookup)
    x[hit] <- lookup[x[hit]]
    x
}
