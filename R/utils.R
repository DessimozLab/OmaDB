source("R/config.R")

depth <- function(list) ifelse(is.list(list), 1L + max(sapply(list, depth)), 0L)

check_response <- function(url){
	resp=httr::GET(url)
	if (httr::http_error(resp)) {
    stop(
      sprintf(
        "OMA API request failed [%s]\n%s", 
        httr::status_code(resp),
        substring(httr::fromJSON(httr::content(resp, "text"))$detail,2) #due to utf-8 encoding
      ),
      call. = FALSE
    )
	}
}

urlGenerator <- function(type,id=NULL,detail=NULL,query_param1=NULL,query_param1_value=NULL,query_param2=NULL,query_param2_value=NULL){
	type = tolower(type)
	
	if(!any(endpoints==type)){
		stop("You must provide a valid endpoint.")
	}

	url_prefix = paste0(API_URL,"/",type,"/")
	if(!is.null(id)){ id=paste0(id,"/")}
	if(!is.null(detail)){ detail=paste0(detail,"/")}
	if(!is.null(query_param1_value)){ query_param1=paste0("?",query_param1,"=",URLencode(query_param1_value))}
	if(!is.null(query_param2_value)){ query_param2=paste0("&",query_param2,"=",URLencode(query_param2_value))}

	final_url= paste0(url_prefix,id,detail,query_param1,query_param2)
	check_response(final_url)
	
	return(final_url)

}

simpleRequest <- function (url){
	response = httr::GET(url)
	content_list = httr::content(response, as = "parsed")
	column_names = names(content_list)

	objectFactory(column_names,content_list)
}


objectFactory <- function(column_names,content_list,type) { 

	list_of_variables = list()
	
	for(name in column_names){
 		list_of_variables[[name]] = content_list[[name]]
 	}

 	value <- list_of_variables
 
 	attr(value, "class") <- "apiobject"
 	value

}

requestFactory <- function (url,type,pagination=TRUE) {

	response = httr::GET(url)
	content_list = httr::content(response, as = "parsed")
	column_names = names(content_list)

	if("count" %in% column_names){
		resolvePagination(url,pagination)
	}

	if(is.null(column_names)){
		formatData(content_list)
	}

	else {
		objectFactory(column_names,content_list,type)
	}
}

resolvePagination <- function(url,pagination) {
	 
	 if(pagination){
		simpleRequest(url)
	 }

	 else{
	 	getAllPages(url)
	 }	
}

get_all_pages <- function(url){
	list_of_things <- list()
	list_of_things[[1]] <- jsonlite::fromJSON(url)$results
	count = jsonlite::fromJSON(url)$count
	pages = count/100
	if (pages%%1 != 0){
		if(pages%%1 >= 0.5){
			pages= round(pages)
		}
		else{
			pages=round(pages+1,digits=0)
		}
	}
	for(i in 2:pages){
		new_url = paste0(url,"?page=",i)
		list_of_things[[i]] <- jsonlite::fromJSON(new_url)$results
	}
	return(list_of_things)
}

