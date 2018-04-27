#' Get the HOG Data Function
#' 
#' The function to obtain the information available for a Hierarchical orthologous group entry in the datase.
#'
#' @param id an identifier for the entry to be returned - either its id or one of its protein members
#' @param level a specific level for the HOG to be restricted to - set to the root level by default. A taxonomic level can be identified by its full capitalised name e.g. "Fungi" or "Saccharomycetaceae". 
#' @param members boolean that when set to TRUE returns a dataframe containg the protein members at a given hog and/or level
#' @return an object containing the JSON keys as attributes
#' @export
#' @examples
#' getHOG(id = "YEAST590")
#' getHOG(id = "YEAST590",level="Saccharomycetaceae", members=TRUE)



getHOG <- function(id, level = NULL, members = FALSE) {
    
    
    if (class(members) != "logical") {
        stop("Members parameter is of type boolean.")
    }
    
    
    if (missing(id)) {
        stop("You must provide a valid HOG ID.")
    }
    
    if (!is.null(level) && class(level) == "integer") {
            stop("You must provide a valid identifier for a taxonomic level - 
                it can only be idenitifed by its full capitalised name.")
        }
    
    if (!is.null(level)) {
        if (members == FALSE) {
            url = urlGenerator(type = "hog", id = id, query_param1 = "level", 
                query_param1_value = level)
        } else {
            url = urlGenerator(type = "hog", id = id, detail = "members", 
                query_param1 = "level", query_param1_value = level)
        }
    } else {
        if (members == FALSE) {
            url = urlGenerator(type = "hog", id = id)
        } else {
            url = urlGenerator(type = "hog", id = id, detail = "members")
        }
    }
    
    return(requestFactory(url))
}









