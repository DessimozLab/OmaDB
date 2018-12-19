
#' Get the Tree Object
#' 
#' A convenience function to obtain a tree object from newick tree,
#' essentially wraps read.tree from the ape package.
#'
#' @param newick The newick tree to be instantiated.
#' @return a tree object
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
