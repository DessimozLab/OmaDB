#' Get the Tree Object
#' 
#' The function to obtain the tree object from newick stored in a text file.
#'
#' @param newick the txt file containing the newick of interest.
#' @return an tree object
#' @export
#' @examples
#' taxonomy = getTaxonomy(root="Alveolata")
#' write(taxonomy$newick, "newick.txt")
#' getTree("newick.txt")


getTree <- function(newick){
	return(ape::read.tree(newick))
}

