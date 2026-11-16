.hogAttributes <- c("members", "gene_ontology", "similar_profile_hogs")

#' Retrieve a HOG from the OMA Browser
#'
#' The function retrieves a specific Hierarchical Orthologous Group (HOG)
#' from the OMA Browser database. A HOG is a set of genes that have all
#' descended from a single ancestral gene at a specific taxonomic level.
#'
#' A HOG can be identified by its member proteins and a taxonomic level, or
#' a HOG ID. As a taxonomic level, you can use either 'root' to retrieve the
#' HOG at its deepest level, or the name of NCBI taxonomy level, or leave it
#' out in which case the deepest level that doesn't include a duplication
#' node is used.
#'
#' The optional argument attribute can be used to directly retrieve
#' \describe{
#'     \item{members}{the member proteins of the HOG at the given level}
#'     \item{gene_ontology}{the Gene Ontology annotations of the HOG, i.e. of
#'     the ancestral gene at the given level, with a score}
#'     \item{similar_profile_hogs}{HOGs with a similar phylogenetic profile,
#'     i.e. similar pattern of gene losses and duplications}
#' }
#'
#' @param id an identifier for the HOG to be returned - either its HOG ID or
#'     a protein id.
#' @param level a specific level for the HOG to be restricted to. level can
#'     either be 'root', or the name of a taxonomic level that is part of the
#'     HOG, e.g. 'Fungi'. By default it will retrieve the deepest level of the
#'     most specific subhog for the given ID.
#' @param members logical, if \code{TRUE} returns a data.frame containing
#'     the member proteins of the HOG at the given level. This is equivalent
#'     to \code{attribute = "members"}.
#' @param attribute an extra attribute of the HOG to be returned, see
#'     details.
#' @param max_results for \code{attribute = "similar_profile_hogs"}, the
#'     maximum number of similar HOGs to return.
#' @return Without \code{attribute}, a data.frame with one row per HOG
#'     matching the query (a HOG ID can refer to several sub-HOGs at a given
#'     level), each with information on the HOG at the given level. With
#'     \code{attribute} (or \code{members = TRUE}), a data.frame with the
#'     requested data.
#' @export
#' @examples
#' getHOG(id = "YEAST00001")
#' getHOG(id = "YEAST00001", level = "root")[, c("hog_id", "level")]
#' members <- getHOG(id = "HOG:F0411479", level = "Saccharomycetaceae",
#'     members = TRUE)
#' head(members$omaid)
#' getHOG(id = "HOG:F0411479", attribute = "gene_ontology")
getHOG <- function(id, level = NULL, members = FALSE, attribute = NULL,
        max_results = NULL) {
    if (missing(id)) {
        stop("You must provide a valid HOG ID.")
    }
    .checkScalarId(id)
    if (!is.logical(members) || length(members) != 1L || is.na(members)) {
        stop("members must be TRUE or FALSE.")
    }
    if (!is.null(level) && !is.character(level)) {
        stop("level must be the name of a taxonomic level or 'root'.")
    }
    if (!is.null(attribute)) {
        attribute <- match.arg(attribute, .hogAttributes)
    } else if (members) {
        attribute <- "members"
    }
    query <- list(level = level)
    if (identical(attribute, "similar_profile_hogs")) {
        query <- list(max_results = max_results)
    }
    res <- .omaGet(.apiUrl("hog", id, attribute, query = query))
    if (identical(attribute, "members")) {
        return(.toDataFrame(res$members))
    }
    if (identical(attribute, "similar_profile_hogs")) {
        return(.toR(res$similar_profile_hogs))
    }
    .toR(res)
}
