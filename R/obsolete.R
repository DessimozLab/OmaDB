#' Defunct functions in package \pkg{OmaDB}
#'
#' These functions were deprecated in OmaDB 1.99.0 and are now defunct.
#' They will be removed in the next release. Use the replacement indicated
#' below:
#' \itemize{
#'     \item \code{getXref}: use \code{\link{searchProtein}}
#'     \item \code{getAnnotation}: use \code{\link{annotateSequence}}
#'     \item \code{getGenomeAlignment}: use \code{\link{getGenomePairs}}
#'     \item \code{getData}: use \code{\link{getProtein}},
#'     \code{\link{getGenome}} or \code{\link{getOMAGroup}}
#' }
#'
#' @param ... ignored
#' @return These functions do not return, they raise an error.
#' @name OmaDB-defunct
#' @aliases getXref getAnnotation getGenomeAlignment getData
#' @keywords internal
NULL

#' @rdname OmaDB-defunct
#' @export
getXref <- function(...) {
    .Defunct("searchProtein", package = "OmaDB")
}

#' @rdname OmaDB-defunct
#' @export
getAnnotation <- function(...) {
    .Defunct("annotateSequence", package = "OmaDB")
}

#' @rdname OmaDB-defunct
#' @export
getGenomeAlignment <- function(...) {
    .Defunct("getGenomePairs", package = "OmaDB")
}

#' @rdname OmaDB-defunct
#' @export
getData <- function(...) {
    .Defunct(msg = paste("'getData' is defunct. Use 'getProtein',",
        "'getGenome' or 'getOMAGroup' instead."), package = "OmaDB")
}
