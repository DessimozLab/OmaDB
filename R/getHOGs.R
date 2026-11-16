#' Retrieve all HOGs at a taxonomic level
#'
#' The function retrieves all Hierarchical Orthologous Groups (HOGs) that
#' are defined at a given taxonomic level. Each HOG corresponds to an
#' ancestral gene of the ancestral genome at that level, so this function
#' can be used to retrieve the gene content of an ancestral genome.
#'
#' If \code{compare_with} is set to a more ancestral taxonomic level, each
#' HOG is annotated with the evolutionary event that happened to the
#' corresponding ancestral gene between the two levels: whether it was
#' retained, duplicated, gained (i.e. it originated between the levels) or
#' lost.
#'
#' @param level the taxonomic level of interest (scientific name or NCBI
#'     taxon id).
#' @param compare_with optional, a more ancestral taxonomic level to compare
#'     with.
#' @return a data.frame with one row per HOG. The column \code{members_url}
#'     can be used with \code{\link{resolveURL}} to retrieve the members of a
#'     HOG. If \code{compare_with} is given, the column \code{event} contains
#'     the evolutionary event.
#' @export
#' @examples
#' hogs <- getHOGs(level = "Saccharomycetaceae")
#' nrow(hogs)
#' head(hogs[, c("hog_id", "completeness_score", "nr_genes")])
getHOGs <- function(level, compare_with = NULL) {
    if (missing(level)) {
        stop("You must provide a taxonomic level.")
    }
    .checkScalarId(level)
    if (!is.null(compare_with)) {
        .checkScalarId(compare_with)
    }
    url <- .apiUrl("hog",
        query = list(level = level, compare_with = compare_with))
    .toR(.omaGetList(url))
}
