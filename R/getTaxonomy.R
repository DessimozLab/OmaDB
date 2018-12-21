#' Get the Taxonomic tree function
#'
#' The function to obtain the taxonomic tree from the database in the newick format that can be plugged into phylo.io for visualisation.
#'
#' @param root optional parameter, the root of the node of interest
#' @param members optional parameter, list of member ncbi taxon or UniProt IDs that should be included in the induced taxonomy.
#' @param newick optional parameter, boolean default set to TRUE
#' @return an object containing the JSON keys as attributes
#' @export
#' @examples
#' getTaxonomy()
#' getTaxonomy(members='YEAST,ASHGO')
#' getTaxonomy(root='Alveolata')

getTaxonomy <- function(root = NULL, members = NULL, newick = TRUE) {
    if (!is.logical(newick)) {
        stop("newick parameter must be a logical value")
    }
    format <- if (newick)
        "newick" else NULL
    url <- urlGenerator(endpoint = "taxonomy", id = root, type = format, members = members)
    return(requestFactory(url))
}
