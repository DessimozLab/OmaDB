#' Retrieve a HOG in OrthoXML format
#'
#' The function retrieves a Hierarchical Orthologous Group (HOG) in the
#' OrthoXML format (\url{https://orthoxml.org}), which contains the full
#' nested structure of the HOG including all duplication and speciation
#' events. The result can for example be analysed with dedicated tools for
#' HOGs such as pyHam.
#'
#' @param id the root HOG ID (e.g. "HOG:F0411479") or the id of a member
#'     protein, in which case the root HOG of the protein is returned.
#' @param augmented logical, whether to return an augmented version of the
#'     OrthoXML that contains more information (e.g. the taxonomic level of
#'     each group), but which does not strictly follow the OrthoXML schema.
#' @return a character string with the OrthoXML document.
#' @export
#' @examples
#' xml <- getOrthoXML("HOG:F0411479")
#' substr(xml, 1, 200)
getOrthoXML <- function(id, augmented = FALSE) {
    if (missing(id)) {
        stop("You must provide a valid HOG ID.")
    }
    .checkScalarId(id)
    if (!startsWith(as.character(id), "HOG:")) {
        id <- getHOG(id, level = "root")$hog_id[1]
    }
    base <- sub("/api/$", "/oma/", .omadb$apiUrl)
    url <- paste0(base, "hog/", utils::URLencode(id, reserved = TRUE),
        "/orthoxml/", if (isTRUE(augmented)) "augmented/")
    .fetch(url, accept = "application/xml, text/plain")$text
}
