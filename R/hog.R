#' An example HOG object. 
#'
#' An object containing information for the HOG:0273533.1b.
#'
#' @format An S3 object with 8 variables:
#' \describe{
#'   \item{hog_id}{hog identifier}
#'   \item{level}{the taxonomic level of this hog}
#'   \item{levels_url}{url pointer to the hog information at a given level}
#'   \item{members_url}{url pointer to the list of gene members for this hog}
#'   \item{alternative_members}{a dataframe object containing the rest of the taxonomic levels in this hog}
#'   \item{roothog_id}{the root taxonomic level of this hog}
#'   \item{parent_hogs}{a dataframe containing information on the parent hogs to the current hogs}
#'   \item{children_hogs}{a dataframe containing information on the children hogs to the current hogs}
#'   ...
#' }
#' @source \url{https://omabrowser.org/api/hog/HOG:0273533.1b/}
"hog"