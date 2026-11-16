#' Get the API and database version
#'
#' The function obtains the version of the OMA Browser database release
#' and of the REST API that the package is using.
#'
#' @return an \code{\link{omadb_obj}} with the attributes
#'     \code{oma_version} and \code{api_version}
#' @export
#' @examples
#' getVersion()
getVersion <- function() {
    .toR(.omaGet(.apiUrl("version")))
}
