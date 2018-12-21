#' Get the API and database version function
#'
#' The function to obtain the API and database version that the package is using.
#'
#' @return S3 object
#' @export
#' @examples
#' getVersion()


getVersion <- function() {
    url <- urlGenerator(endpoint = "version")
    return(requestFactory(url))
}
