#' Get the BioString AAString Object
#' 
#' The function to obtain the information available for a single entry in the datase.
#'
#' @param AAString the amino acid string of interest
#' @return an AAString object
#' @export
#' @examples
#' getAAString(AAString="NPROK")


getAAString <- function(AAString){
	if(class(AAString)=="character"){
		return(Biostrings::AAString(AAString))
	}
	else{
		stop("The amino acid string must be a string of characters.")
	}
	
}



#' Get the BioString DNAString Object
#' 
#' The function to obtain the information available for a single entry in the datase.
#'
#' @param DNAString the DNA string of interest
#' @return an DNAString object
#' @export
#' @examples
#' getDNAString(DNAString="ATGC")

getDNAString <- function(DNAString){
	if(class(DNAString)=="character"){
		return(Biostrings::DNAString(DNAString))
	}
	else{
		stop("The DNA string must be a string of characters.")
	}
	
}

