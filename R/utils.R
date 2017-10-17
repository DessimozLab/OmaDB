
API_URL = "https://omabrowser.org/api"

depth <- function(list) ifelse(is.list(list), 1L + max(sapply(list, depth)), 0L)


#' @importFrom utils URLencode
#' @import httr
#' @import plyr
#' @importFrom  Biostrings AAString
#' @importFrom Biostrings DNAString 
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom IRanges IRanges

urlGenerator <- function(type = NULL, id = NULL, detail = NULL, query_param1 = NULL, 
    query_param1_value = NULL, query_param2 = NULL, query_param2_value = NULL, 
    query_param3 = NULL, query_param3_value = NULL) {
    
    if (!is.null(type)) {
        type = tolower(type)
    }
    
    
    url_prefix = paste0(API_URL, "/", type, "/")
    if (!is.null(id)) {
        id = paste0(id, "/")
    }
    if (!is.null(detail)) {
        detail = paste0(detail, "/")
    }
    
    if (!is.null(query_param1_value)) {
        query_param1 = paste0("?", query_param1, "=", utils::URLencode(query_param1_value))
    }
    if (!is.null(query_param2_value)) {
        query_param2 = paste0("&", query_param2, "=", utils::URLencode(query_param2_value))
    }
    if (!is.null(query_param3_value)) {
        query_param3 = paste0("&", query_param3, "=", utils::URLencode(query_param3_value))
    }
    
    final_url = paste0(url_prefix, id, detail, query_param1, query_param2, 
        query_param3)
    
    return(final_url)
    
}

simpleRequest <- function (url){


	out <- tryCatch(
	{		
		response = httr::GET(url)

		content_list = httr::content(response, as = "parsed")
		column_names = names(content_list)

		objectFactory(column_names,content_list)
	}, 
	error = function(cond) {
            message(paste("THE OMA REST API request failed:", url))
            message("Here's the original error message:")
            
			response_message = httr::http_status(response)$message
			
            message(response_message)
            
            return(NA)
        },

    warning = function(cond) {
            message(paste("URL caused a warning:", url))
            message("Here's the original warning message:")
            
            response_message = httr::http_status(response)$message

            message(response_message)
 
            return(NULL)
        }
   
   	)

	return(out)
	
}


objectFactory <- function(column_names, content_list) {
    
    list_of_variables = list()
    
    for (name in column_names) {
        
        content = content_list[[name]]
        
        if (class(content) == "list" && length(content)!=0) {
            if (is.null(names(content))) {
                  content = formatData(content)
            }
        }
        
        if (name == "chromosomes") {
            content = GenomicRanges::makeGRangesFromDataFrame(content, 
                start.field = "entry_ranges.1", end.field = "entry_ranges.2", 
                seqnames.field = "id", ignore.strand = TRUE)
        }
        
        if (name == "locus") {
            content = GenomicRanges::GRanges(seqnames = content_list[["omaid"]], 
                ranges = IRanges::IRanges(content$start, content$end), 
                strand = content$strand)
        }
        
        if (name == "sequence") {
            content = Biostrings::AAString(content)
        }
        
        if (name == "cdna" && !(grepl("X", content))) {
            content = Biostrings::DNAString(content)
        }
        
        list_of_variables[[name]] = content
    }
    
    value <- list_of_variables
    
    return(value)
    
}

requestFactory <- function (url) {

	

	out<- tryCatch(
	{	
		response = httr::GET(url)

		content_list = httr::content(response, as = "parsed")
		column_names = names(content_list)

		if(is.null(column_names)){
			if(length(content_list)==1){
				column_names = names(content_list[[1]])

				objectFactory(column_names, content_list = content_list[[1]])
				}
			
			else{

				formatData(content_list)
			}
			
		}

		else {

			objectFactory(column_names,content_list)
		}
	},
		
		error= function(cond) {
            message(paste("THE OMA REST API request failed: the endpoint does not exist:", url))
            message("Here's the original error message:")
            
			response_message = httr::http_status(response)$message
			
            message(response_message)
            
            return(NA)
        },

    	warning=function(cond) {
            message(paste("URL caused a warning:", url))
            message("Here's the original warning message:")
            
            response_message = httr::http_status(response)$message

            message(response_message)
 
            return(NULL)
        }	

		)
	
	return(out)
}

formatData <- function(data) {
    if (class(data) == "list" && length(data) != 0) {
        
        if ("entry_1" %in% names(data[[1]])) {
            for (i in seq_along(data)) {
                data[[i]][[1]][[7]] = rbind(data[[i]][[1]][[7]])
                data[[i]][[2]][[7]] = rbind(data[[i]][[2]][[7]])
            }
            
        }
        
        if ("entry_ranges" %in% names(data[[1]])) {
            
            for (i in seq_along(data)) {
                data[[i]][[2]][[1]] = rbind(data[[i]][[2]][[1]])
                
            }
            
        }
        
        
        dfs <- lapply(data, data.frame, stringsAsFactors = FALSE)
        data = plyr::rbind.fill(dfs)
        
        return(data)
        
    }
    
}








