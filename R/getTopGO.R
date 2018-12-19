#' Get the topGO Object function
#' 
#' The function to create a topGO object containing the GO annotations for the given protein list. 
#'
#' @param annotations list of GO annoatations obtained from the formatTopGO()
#' @param foregroundGenes List of identifiers for the foreground genes
#' @param format Format for the data to be returned in - either 'GO2geneID' or 'geneID2GO'
#' @param ontology The ontology for which the enrichment should be done. This 
#'        parameter is passed directly to the topGOdata constructor.
#' @return topGO object
#' @export
#' @import topGO
#' @importFrom topGO annFUN.gene2GO
#' @importFrom methods new
#' @examples 
#' geneList = list(getData(type="protein",id="YEAST58"),getData(type="protein",id="YEAST00059"))
#' annotations = formatTopGO(geneList,format="geneID2GO")
#' library(topGO)
#' getTopGO(annotations, foregroundGenes = list("YEAST00058"), format = "geneID2GO")

getTopGO <- function(annotations, format, foregroundGenes, ontology){
	
	if(missing(annotations) || class(annotations) != "list"){
		stop("You must provide a valid list of annotations.")
	}

	if(missing(format) || !(format %in% list("geneID2GO","GO2geneID"))){
		stop("You must provide a valid annotations format.")
	}

	if(missing(foregroundGenes)){
		stop("You must provide a valid list of genes of interest.")
	}

	if(class(foregroundGenes)=="data.frame"){
		foregroundGenes = foregroundGenes[['omaid']]
	}

	geneNames <- names(annotations)
	geneList <- factor(as.integer(geneNames %in% foregroundGenes))
	names(geneList) <- geneNames

	if(format=="geneID2GO"){
		GOdata = new("topGOdata", ontology = ontology, allGenes = geneList,
              annot = topGO::annFUN.gene2GO, gene2GO = annotations)
	}
	if(format=="GO2geneID"){
		GOdata = new("topGOdata", ontology = ontology, allGenes = geneList,
              annot = topGO::annFUN.GO2genes, GO2gene = annotations)
	}
	return(GOdata)
}



