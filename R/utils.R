
API = "https://omabrowser.org/api/"


depth <- function(list) ifelse(is.list(list), 1L + max(sapply(list, depth)), 0L)


#' @importFrom utils URLencode
#' @import httr
#' @import plyr
#' @import ape
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
    
    
    url_prefix = paste0(API, type, "/")
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


simpleRequest <- function (url,body = NULL){
    if(is.null(body)){
        response = load_url(url)
    }
    else{
        response = load_url(url,body)
    }

    if (!is.null(response)){
        content_list = httr::content(response, as = "parsed")
        column_names = names(content_list)
        return( objectFactory(column_names,content_list) )
    }
    return(NULL)
}

load_url <- function(url, body = NULL){


    count = 0
    while(count<3){
        response <- tryCatch(
            if(!is.null(body)){
                 httr::POST(url, body = body, encode = "raw", accept('application/json'),content_type('application/json'))

            }
            else {
                httr::GET(url)
            }
            , 
            error = function(cond){
                Sys.sleep(0.5)
                return(NULL)
            },
            warning = function(cond){
                message(cond)
                return(NULL)
            }
        )

        if (is.null(response)){
            count = count + 1
            next
        }

        if (response$status_code < 500){
            break
        }
        message(sprintf("Request failed with server error %d. retry in %d secs", response$status_code, 2**count))
        Sys.sleep(2**count)
        count = count + 1
    }

    if (is.null(response)){
        message(sprintf("Cannot resolve '%s'. Please connect to the internet", url))
        return(NULL)
    }
    
    tryCatch(
    {
        httr::stop_for_status(response)
        # request worked out 
        return(response)
    }, 
    error = function(cond) {
        message("THE OMA REST API request failed:", url)
        message("Here's the original error message:")
        message(cond)
        
        return(NULL)
    },
    warning = function(cond) {
        message("URL caused a warning:", url)
        message("Here's the original warning message:")
        message(cond)
    
        return(NULL)
    }
    )
}


objectFactory <- function(column_names, content_list) {
    

    list_of_variables = lapply(column_names, FUN = function(name) {

            content = content_list[[name]]

            if(is.null(content)){
                content == " "
            }

            if (class(content) == "list" && length(content)!=0 && name!="locus"  && name!="chromosomes") {
                if (is.null(names(content))) {
                    formatData(content)  
                }
            }

            else if (name == "chromosomes") {
                content = formatData(content)
                GenomicRanges::makeGRangesFromDataFrame(content, 
                    start.field = "entry_ranges.1", end.field = "entry_ranges.2", 
                    seqnames.field = "id", ignore.strand = TRUE)
            }
            
            else if (name == "locus") {
                GenomicRanges::GRanges(seqnames = content_list[["omaid"]], 
                    ranges = IRanges::IRanges(content$start, content$end), 
                    strand = content$strand)
            }
            
            else if (name == "sequence") {
                Biostrings::AAString(content)
            }
            
            else if (name == "cdna") {
                content = gsub('X','N', content)
                Biostrings::DNAString(content)
            }
            
            else{
                content
            }

        })


    names(list_of_variables) = column_names
    
    value <- list_of_variables

    class(value) <- 'omadb_obj'
    
    return(value)
    
}

largeRequestFactory <- function(url, n) {


    n_requests = round(as.numeric(n)/10000)

    url_list = list()

    for(i in seq_along(1:n_requests)){

        url_list[[i]] = paste0(url,'?per_page=10000&page=',i)

    }

    response_list = lapply(url_list, FUN = function(x) { load_url(x) } )

    json_list = lapply(response_list, FUN = function(x) { httr::content(x,as = 'parsed') } )

    content_list = do.call("c",json_list)

    column_names = names(content_list)

    if(is.null(column_names)){
        if(length(content_list)==1){
            column_names = names(content_list[[1]])
            content_list = content_list[[1]]

            return(objectFactory(column_names, content_list))
         
        }  

        else if (length(content_list)!=1 && length(content_list)!=0){

            return(formatData(content_list))

            }
   
        else if (length(content_list)==0){

            return(" ")

        }
            
    } 

    else {

        return( objectFactory(column_names,content_list))

    }

return(NULL)
}

requestFactory <- function (url,body=NULL,per_page=NULL) {

    
    if(is.null(body)){
        response = load_url(url)
    }
    else{
        response = load_url(url,body)
    }

    if(!is.null(per_page)){
        if(per_page=='all'){
            n_items = headers(response)[['x-total-count']]
            
            if(as.numeric(n_items)>10000){
                largeRequestFactory(url,n = n_items)
            }

            new_url = paste0(url,"?per_page=",as.character(n_items))
        }
        else{

            if(per_page>10000){
                largeRequestFactory(url,n = per_page)
            }

            new_url = paste0(url,"?per_page=",as.character(per_page))
            
        }

        response = load_url(new_url)
    }

    if (!is.null(response)){
        content_list = httr::content(response, as = "parsed")
        column_names = names(content_list)

        if(is.null(column_names)){
            if(length(content_list)==1){
                column_names = names(content_list[[1]])
                content_list = content_list[[1]]

                return(objectFactory(column_names, content_list))
             
            }  

            else if (length(content_list)!=1 && length(content_list)!=0){
                return(formatData(content_list))
                }
       
            else if (length(content_list)==0){
                return(" ")
            }
                
        } else {
            return( objectFactory(column_names,content_list))
        }
    }
    return(NULL)
}

formatData <- function(data) {
        ## whole genome alignment
        if ("entry_1" %in% names(data[[1]])) {
            for (i in seq_along(data)) {
                data[[i]][[1]][[7]] = rbind(data[[i]][[1]][[7]])
                data[[i]][[2]][[7]] = rbind(data[[i]][[2]][[7]])
            }
            
        }
        ##hogs
        else if("alternative_levels" %in% names(data[[1]])){
            for (i in seq_along(data)) {
                data[[i]][['alternative_levels']] = NULL
            }        
        }
        ## flatten loci
        else if ("entry_ranges" %in% names(data[[1]])) {
            
            for (i in seq_along(data)) {
                data[[i]][['entry_ranges.1']] = data[[i]]['entry_ranges'][[1]][[1]][[1]]
                data[[i]][['entry_ranges.2']] = data[[i]]['entry_ranges'][[1]][[1]][[2]]
                data[[i]][['entry_ranges']] = NULL
            }
            
        }
        
        else if("sequence" %in% names(data[[1]])){
            for (i in seq_along(data)) {
                data[[i]] = objectFactory(names(data[[i]]),data[[i]])
                
            }
            return(data)
        }

     
        dfs <- lapply(data, data.frame, stringsAsFactors = FALSE)
        data = plyr::rbind.fill(dfs)
        
        return(data)
    
}

#' Resolve URLs automatically when accessed
#' 
#' The function to obtain further information from a given url. 
#'
#' @param x object
#' @param name attribute
#' @return API response behind the URL
#' @export

'$.omadb_obj' <- function(x,name) {

    if(grepl('https://',x[[name]])){

        value <- resolveURL(x[[name]])
        obj_name = deparse(substitute(x))
        set_new_val(x, name) = value
        assign(obj_name, x, envir = .GlobalEnv)

    }

    return(x[[name]])

}

"set_new_val<-" = function(x, name, value) {

    x[[name]] = value
    x
}












