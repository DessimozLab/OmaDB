#' Get the Taxonomic tree function
#' 
#' The function to obtain the taxonomic tree from the database in the newick format that can be plugged into phylo.io for visualisation. 
#'
#' @param members optional parameter, list of member ncbi taxon or UniProt IDs that should be included in the induced taxonomy. 
#' @param root optional parameter, the root of the node of interest
#' @param newick optional parameter, boolean default set to TRUE
#' @return an object containing the JSON keys as attributes
#' @export
#' @examples
#' getTaxonomy()
#' getTaxonomy(members="YEAST,ASHGO")
#' getTaxonomy(root="Alveolata")

getTaxonomy <- function(root=NULL,members,newick=TRUE) {
	if(missing(members)){
		if(newick==FALSE){
			url = urlGenerator(type="taxonomy",id=root)

		}
		else{
			url = urlGenerator(type="taxonomy",id=root,query_param1="type",query_param1_value="newick")

		}
	}
	else{
		if(newick==FALSE){
			url = urlGenerator(type="taxonomy",id=root,query_param1="members",query_param1_value=members)
		}
		else{
			url = urlGenerator(type="taxonomy",id=root,query_param1="members",query_param1_value=members,
				query_param2="type",query_param2_value="newick")
		}
						
	}
	return(requestFactory(url))
}