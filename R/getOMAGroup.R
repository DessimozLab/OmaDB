#' Retrieve an OMA Group from the OMA Browser
#'
#' This function obtains an OMA Group from the OMA Browser database. An OMA Group
#' is defined to be a clique of proteins that are all orthologous to each other,
#' i.e. they are all related through speciation events only. An OMA Group can thus
#' by definition not contain any inparalogs. It is a very stringent orthology
#' grouping approach.
#' OMA Groups are mostly useful to infer phylogenetic species tree where they
#' can be used as marker genes.
#'
#' Retrieving an OMA Group can be done using a group nr as id, its fingerprint
#' (a 7mer AA sequence which is unique to proteins in that group), a member
#' protein id or any sequence pattern that is unique to the group.
#'
#' @param id An identifier for the group. See above for possible types of IDs.
#' @param attribute an extra attribute to be returned (close_groups)
#' @return an object containing the JSON keys as attributes or a dataframe
#' @export
#' @examples
#' getOMAGroup(id='58')
#' getOMAGroup(id='P12345')
#' getOMAGroup(id='NNRRGRI')
#' getOMAGroup(id='58', attribute='close_groups')

getOMAGroup <- function(id, attribute = NULL) {

    if (!is.null(attribute) && !(attribute %in% c("close_groups"))) {
        stop("You must provide a valid attribute.")
    }
    url <- urlGenerator(endpoint = "group", id = id, detail = attribute)
    return(requestFactory(url))
}
