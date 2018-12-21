#' Retrieve a HOG from the OMA Browser
#'
#' The function retrieves a specific Hierarchical Orthologous Group (HOG) from
#' the OMA Browser database. A HOG is a set of genes that have all decendet from
#' a single ancestral gene at a specific taxonomic level.
#'
#' A HOG can be identified by its member proteins and a taxonomic level, or
#' a HOG ID. As a taxonomic level, you can use either 'root' to retrieve the
#' HOG at its deepest level, or the name of NCBI taxonomy level, or leave it
#' out in which case the deepest level that doesn't include a duplication node
#' is used.
#'
#' The function either returns a single hog object or a list of hog
#' objects. The later happens if the HOG ID you provide has already split into
#' several sub-hogs at the level you indicate.
#'
#' @param id an identifier for the HOG to be returned - either its HOG ID or a protein id.
#' @param level a specific level for the HOG to be restricted to. level can either be 'root', or the name of a taxonomic level that is part of the HOG, e.g. 'Fungi'. By default it will retrieve the depest level of the most specific subhog for the given ID.
#' @param members boolean that when set to TRUE returns a dataframe containg the protein members at a given hog level
#' @return an object containing HOG attributes, or a list of those
#' @export
#' @examples
#' getHOG(id = 'YEAST590')
#' getHOG(id = 'YEAST590', level='root')
#' getHOG(id = 'YEAST590', level='Saccharomycetaceae', members=TRUE)

getHOG <- function(id, level = NULL, members = FALSE) {

    if (class(members) != "logical") {
        stop("Members parameter is of type boolean.")
    }

    if (missing(id)) {
        stop("You must provide a valid HOG ID.")
    }

    if (!is.null(level) && class(level) == "integer") {
        stop("You must provide a valid identifier for a taxonomic level -
                it can only be idenitifed by its full capitalised name.")
    }

    memb <- if (members)
        "members" else NULL
    url <- urlGenerator(endpoint = "hog", id = id, detail = memb, level = level)
    return(requestFactory(url))
}









