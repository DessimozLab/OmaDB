#' Get further information for a dataframe of members
#' 
#' The function to obtain further information from a dataframe containing a list of members. 
#'
#' @param df the dataframe or a list of dataframes containing the genomic range data of interest
#' @param type the type of information to be retrieved
#' @return an list 
#' @export
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @examples
#' sequences = getInfo(df = getData("group","YEAST58")['members'],type='sequences')


getInfo <- function(df,type){

	if(!(type %in% c('genomic_ranges','domains','sequences','ontologies'))){
		stop("You must enter a valid type of information to be retrieved.")
	}

	if(class(df)=="list"){
			df = plyr::rbind.fill(df)
		}

	if(type == 'domains'){

		response = lapply(df$omaid, FUN = function (x) getData(type = "protein", x,'domains'))

	}

	if(type == 'ontologies'){

		object_list = lapply(df$omaid, FUN = function (x) getData(type = "protein", x))

		response = formatTopGO(object_list,format="geneID2GO")
	}
	if(type == 'sequences'){

		object_list = lapply(df$omaid, FUN = function (x) getData(type = "protein", x))

		response = lapply(object_list, FUN = function(x) x$sequence)

		names(response) = df$omaid
	}
	if(type == 'genomic_ranges'){

		df$locus.strand[df$locus.strand == "1"] <- "+"
		df$locus.strand[df$locus.strand == "-1"] <- "-"

		response = GenomicRanges::makeGRangesFromDataFrame(df,start.field = "start", 
			end.field = "end", seqnames.field="omaid")
	}

	return(response)

}



