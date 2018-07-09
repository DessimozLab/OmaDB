
#' Get the Tree Object
#' 
#' The function to obtain the tree object from newick.
#'
#' @param newick the newick of interest.
#' @return an tree object
#' @export
#' @importFrom ape read.tree
#' @examples
#' taxonomy = getTaxonomy(root="Alveolata")
#' getTree(newick=taxonomy$newick)


getTree <- function(newick){
	if(class(newick)=="list"){
		newick = newick$newick
	}
	return(ape::read.tree(text=newick))
}