
#' Format the GO annotations data
#' 
#' The function to create a list of GO annotations that is compatible with topGO from protein objects in roma 
#'
#' @param geneList the list of roma protein objects to be included in the analysis - this is where the GO annotations are extracted from
#' @param format format for the data to be returned in - either 'GO2geneID' or 'geneID2GO'
#' @return either GO2geneID list or geneID2GO depending on the format parameter
#' @export
#' @examples
#' geneList = list(getData(type="protein",id="YEAST58"),getData(type="protein",id="YEAST00059"))
#' annotations = formatTopGO(geneList,format="geneID2GO")



formatTopGO <- function(geneList,format){

	geneID2GO = list()

	for(protein in geneList){
		
		if(startsWith(protein['ontology'][[1]],"https://")){
			protein['ontology'][[1]] = resolveURL(protein['ontology'][[1]])
		}

		protein_annotations = {}

		for (GO_term in protein['ontology'][[1]]['GO_term']){
			protein_annotations <- GO_term 
		}

		geneID2GO[[protein$omaid]] = protein_annotations
	}

	if(format=="geneID2GO"){
			return(geneID2GO)
			}
	if(format=="GO2geneID"){
			return(topGO::inverseList(geneID2GO))
			}
	else{
		stop("Invalid format. Must be either 'GO2geneID' or 'geneID2GO'")
	}

}

#' Get the topGO Object function
#' 
#' The function to create a topGO object containing the GO annotations for the given protein list. 
#'
#' @param annotations list of GO annoatations obtained from the formatTopGO()
#' @param myInterestingGenes list identifiers for the genes of interest
#' @param format format for the data to be returned in - either 'GO2geneID' or 'geneID2GO'
#' @return topGO object
#' @export
#' @examples
#' annotations = formatTopGO(geneList = list(getData(type="protein",id="YEAST58"),getData(type="protein",id="YEAST00059")),format="geneID2GO")
#' library(topGO)
#' getTopGO(annotations=annotations, myInterestingGenes = list("YEAST00058"), format = "geneID2GO")




getTopGO <- function(annotations,format,myInterestingGenes){
	
	geneNames <- names(annotations)

	geneList <- factor(as.integer(geneNames %in% myInterestingGenes))

	names(geneList) <- geneNames

	if(format=="geneID2GO"){
		GOdata = new("topGOdata", ontology = "MF", allGenes = geneList,
              annot = topGO::annFUN.gene2GO, gene2GO = annotations)
	}
	if(format=="GO2geneID"){
		GOdata = new("topGOdata", ontology = "MF", allGenes = geneList,
              annot = topGO::annFUN.GO2gene, GO2gene = annotations)
	}

	return(GOdata)
	
}



