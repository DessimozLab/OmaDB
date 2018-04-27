#' An example OMA group object. 
#'
#' An object containing information for the OMA group number 737636.
#'
#' @format An S3 object with 4 variables:
#' \describe{
#'   \item{group_nr}{group number, not stable across releases}
#'   \item{fingerprint}{fingerprint of the oma group, stable across releases}
#'   \item{related_groups}{url to the endpoint containing the list of oma groups that share some of the orthologs with this oma group}
#'   \item{members}{list of protein members of this oma group}
#'   ...
#' }
#' @source \url{https://omabrowser.org/api/group/YEAST58/}
"group"