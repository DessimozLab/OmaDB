#' Bulk retrieve information for a list of proteins
#' 
#' The function to bulk retrieve information for a list of proteins.
#'
#' @param protein_list list of protein members
#' @return a list of S3 objects
#' @export
#' @examples
#' orthologs = getData(type="protein",id='YEAST58')$orthologs
#' bulkRetrieve(orthologs)

#' @importFrom  jsonlite toJSON


bulkRetrieve <- function(protein_list){
	body = jsonlite::toJSON(list(ids=protein_list,auto_unbox=TRUE))
	url = "api/protein/bulk_retrieve/"
	return(requestFactory(url = url, body = body))
}