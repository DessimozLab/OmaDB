#' Get the Object Attributes
#' 
#' The function to obtain the attributes and their data types for the object created.
#'
#' @param obj the object of interest
#' @return an list of object attributes and their data classes
#' @export
#' @examples
#' attributes = getObjectAttributes(getData("group","YEAST58"))


getObjectAttributes <- function (obj){
	for(name in attributes(obj)$names){

 		class_name = class(obj[[name]])

 		if(class_name=="character" && startsWith(obj[[name]],'https://')){
 			class_name = "URL"

 		}
 		
 		print(paste(name, ":", class_name))
 		
 	}
}

#' Get the value for the Object Attribute
#' 
#' The function to obtain the value for an object attribute.
#'
#' @param obj the object of interest
#' @param attribute the attribute of interest
#' @return an value for a given object attribute
#' @export
#' @examples
#' members = getAttribute(getData("group","YEAST58"),'members')


getAttribute <- function (obj,attribute){
	return(obj[[attribute]])
}