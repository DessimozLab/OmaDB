#' Get the GRanges object
#' 
#' The function to obtain the GRanges object from a dataframe containing a list of members. 
#'
#' @param df the dataframe or a list of dataframes containing the genomic range data of interest
#' @return an GRanges object
#' @export
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @examples
#' gRanges = getGRanges(df = getData("group","YEAST58")['members'])



getGRanges <- function(df){

	if(class(df)=="list"){
		df = plyr::rbind.fill(df)
	}

	df$locus.strand[df$locus.strand == "1"] <- "+"
	df$locus.strand[df$locus.strand == "-1"] <- "-"

	grange = GenomicRanges::makeGRangesFromDataFrame(df,start.field = "start", 
		end.field = "end", seqnames.field="omaid")

	return(grange)
}

#' Get the gene ontologies from members
#' 
#' The function to obtain the gene ontologies from a dataframe of members.
#'
#' @param df the dataframe or a list of dataframes containing the member proteins of interest
#' @return a list containing the geneID2GO information 
#' @export
#' @examples
#' ontologies = getOntologies(df = getData("group","YEAST58")['members'])


getOntologies <- function(df){

	if(class(df)=="list"){
		df = plyr::rbind.fill(df)
	}

	object_list = lapply(df$omaid, FUN = function (x) getData(type = "protein", x))

	annotations = formatTopGO(object_list,format="geneID2GO")

	return(annotations)

}

#' Get the sequences from members
#' 
#' The function to obtain the protein seqeunces from a dataframe of members.
#'
#' @param df the dataframe or a list of dataframes containing the member proteins of interest
#' @return a list containing the AAString object for each member gene
#' @export
#' @examples
#' sequences = getSequences(df = getData("group","YEAST58")['members'])

getSequences <- function(df){

	if(class(df)=="list"){
		df = plyr::rbind.fill(df)
	}

	object_list = lapply(df$omaid, FUN = function (x) getData(type = "protein", x))

	sequences = lapply(object_list, FUN = function(x) x$sequence)

	names(sequences) = df$omaid
  	
	return(sequences)

}