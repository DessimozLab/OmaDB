
API_URL = "https://omabrowser.org/api"

depth <- function(list) ifelse(is.list(list), 1L + max(sapply(list, depth)), 0L)

urlGenerator <- function(type,id=NULL,detail=NULL,query_param1=NULL,query_param1_value=NULL,query_param2=NULL,query_param2_value=NULL,query_param3=NULL,query_param3_value=NULL){
	type = tolower(type)
	
	url_prefix = paste0(API_URL,"/",type,"/")
	if(!is.null(id)){ id=paste0(id,"/")}
	if(!is.null(detail)){ detail=paste0(detail,"/")}
	if(!is.null(query_param1_value)){ query_param1=paste0("?",query_param1,"=",utils::URLencode(query_param1_value))}
	if(!is.null(query_param2_value)){ query_param2=paste0("&",query_param2,"=",utils::URLencode(query_param2_value))}
	if(!is.null(query_param3_value)){ query_param3=paste0("&",query_param3,"=",utils::URLencode(query_param3_value))}

	final_url= paste0(url_prefix,id,detail,query_param1,query_param2,query_param3)

	return(final_url)

}

simpleRequest <- function (url){

	response = httr::GET(url)

	out <- tryCatch(
	{		
		content_list = httr::content(response, as = "parsed")
		column_names = names(content_list)

		objectFactory(column_names,content_list)
	}, 
	error= function(cond) {
            message(paste("THE OMA REST API request failed:", url))
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


objectFactory <- function(column_names,content_list) { 

	list_of_variables = list()
	
	for(name in column_names){
 		list_of_variables[[name]] = content_list[[name]]
 	}

 	value <- list_of_variables
 
 	attr(value, "class") <- "OMAObject"
 	value

}

requestFactory <- function (url) {

	response = httr::GET(url)

	out<- tryCatch(
	{
		content_list = httr::content(response, as = "parsed")
		column_names = names(content_list)

		if(is.null(column_names)){
			if(length(content_list)==1){
				column_names = names(content_list[[1]])

				objectFactory(column_names,content_list[[1]])
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
            message(paste("THE OMA REST API request failed:", url))
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


check_response <- function(url){

		response = httr::GET(url)

		if(httr::http_status(response)$category != "Success"){
			message = httr::http_status(response)$message
		
			stop(paste("THE OMA REST API request failed:", url,"\nHere's the original error message:", message))

		}
}


