
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

#' Get the topGO Object function
#' 
#' The function to create a topGO object containing the GO annotations for the given protein list. 
#'
#' @param annotations list of GO annoatations obtained from the formatTopGO()
#' @param myInterestingGenes list of identifiers for the genes of interest or a dataframe containing them
#' @param format format for the data to be returned in - either 'GO2geneID' or 'geneID2GO'
#' @return topGO object
#' @export
#' @import topGO
#' @importFrom topGO annFUN.gene2GO
#' @importFrom methods new
#' @examples 
#' geneList = list(getData(type="protein",id="YEAST58"),getData(type="protein",id="YEAST00059"))
#' annotations = formatTopGO(geneList,format="geneID2GO")
#' library(topGO)
#' getTopGO(annotations, myInterestingGenes = list("YEAST00058"), format = "geneID2GO")




getTopGO <- function(annotations,format,myInterestingGenes){
	
	if(missing(annotations) || class(annotations) != "list"){
		stop("You must provide a valid list of annotations.")
	}

	if(missing(format) || !(format %in% list("geneID2GO","GO2geneID"))){
		stop("You must provide a valid annotations format.")
	}

	if(missing(myInterestingGenes)){
		stop("You must provide a valid list of genes of interest.")
	}

	if(class(myInterestingGenes)=="data.frame"){
		myInterestingGenes = myInterestingGenes[['omaid']]
	}

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



