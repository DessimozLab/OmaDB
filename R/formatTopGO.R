
#' Format the GO annotations data
#' 
#' The function to create a list of GO annotations that is compatible with topGO from protein objects in roma 
#'
#' @param geneList the list of roma protein objects to be included in the analysis - this is where the GO annotations are extracted from
#' @param format format for the data to be returned in - either 'GO2geneID' or 'geneID2GO'
#' @return a list containing the GO2geneID or geneID2GO information
#' @export
#' @examples
#' geneList = list(getData(type="protein",id="YEAST01"),getData(type="protein",id="YEAST03"))
#' annotations = formatTopGO(geneList,format="geneID2GO")



formatTopGO <- function(geneList,format){

	if(!(format %in% list("GO2geneID","geneID2GO"))){
		stop("Invalid format. Must be either 'GO2geneID' or 'geneID2GO'")
	}

	geneID2GO = lapply(geneList, FUN = function(protein) {
		if(startsWith(protein[['ontology']],"https://")){
			annotation = resolveURL(protein[['ontology']])
			if(class(annotation)=="data.frame"){
				unlist(as.list(annotation[["GO_term"]]))
			}
			else{
				annotation
			}		
			}
		else{
			if(class(annotation)=="data.frame"){
				unlist(as.list(protein[['ontology']][['GO_term']]))
			}
			else{
				annotation
			}			
		}
		})

	names(geneID2GO) = lapply(geneList, FUN = function(protein){ protein$omaid })

	if(format=="geneID2GO"){
			return(geneID2GO)
			}
	if(format=="GO2geneID"){
			return(topGO::inverseList(geneID2GO))
			}
}