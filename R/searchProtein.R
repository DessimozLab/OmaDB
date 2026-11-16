#' Search for proteins by a (partial) identifier
#'
#' The function lists all the cross-references to other databases (e.g.
#' UniProtKB, RefSeq, Ensembl) of proteins in the OMA Browser that match a
#' pattern.
#'
#' @param pattern the pattern to query the OMA database with - needs to be
#'     at least 3 characters long. If the pattern matches too many
#'     cross-references, an error is raised.
#' @return a data.frame containing information on the cross references
#'     matching the pattern, including the OMA id of the protein (column
#'     \code{omaid}).
#' @export
#' @examples
#' res <- searchProtein(pattern = "P53_RAT")
#' res[, c("xref", "source", "omaid")]
searchProtein <- function(pattern) {
    if (missing(pattern) || !is.character(pattern) ||
            length(pattern) != 1L || nchar(pattern) < 3L) {
        stop("pattern must be a string with at least 3 characters.")
    }
    res <- .omaGet(.apiUrl("xref", query = list(search = pattern)))
    if (.isScalarList(res) && length(res) > 0L) {
        stop(paste(unlist(res), collapse = " "))
    }
    .toR(res)
}
