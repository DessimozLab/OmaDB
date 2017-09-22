#' Get the Formatted Data Function
#' 
#' The function to convert a named list into a corresponding dataframe. 
#'
#' @param data the list to be formatted
#' @return a dataframe
#' @export
#' @examples
#' formatData(getData(type="Genome",id="YEAST")$chromosomes)

formatData <- function (data) {
	if(class(data)=="list" && length(data)!=0){
	 		if(depth(data)>3){ #chromosomes
	 			for(i in 1:length(data)){
	 				data[[i]][[2]][[1]]=rbind(data[[i]][[2]][[1]])		
	 			}

	 		}

	 		dfs <- lapply(data, data.frame, stringsAsFactors = FALSE)
	 		data = plyr::rbind.fill(dfs)
	 		
	 		return(data)
        
	}

}