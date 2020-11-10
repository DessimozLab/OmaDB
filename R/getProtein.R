#' Retrieve a protein from the OMA Browser
#'
#' This function enables to retrieve information on one or several proteins from
#' the OMA Browser database.
#'
#' In its simplest form the function returns the base data of the query protein.
#' The query protein can be selected with any unique id, for example with a
#' UniProtKB accession (P12345), an OMA id (YEAST00012), or a RefSeq id
#' (NP_001226). To retrieve more than one protein, you should pass a vector of
#' IDs.
#'
#' Non-scalar properties of proteins such as their domains, GO annotations,
#' orthologs or homeologs will get loaded upon accessing them, or if you only
#' need this information you can set the attribute parameter to the property
#' name and retrieve this information directly.
#'
#' @seealso For non-unique non-unique IDs or partial ID lookup, use
#'   \code{\link{searchProtein}} instead.
#'
#' @param id Identifier(s) for the entry or entries to be returned. a character
#'   string if single entry or a vector if multiple.
#' @param attribute Instead of the protein, return the attribute property of the
#'   protein. Attriute needs to be one of 'domains', 'orthologs',
#'   'gene_ontology', 'locus', or 'homoeologs'.
#' @return An object containing the JSON keys as attributes or a dataframe
#'   containing the non-scalar protein property.
#' @export
#' @examples
#' getProtein(id='YEAST00001')
#' getProtein(id='YEAST00001', attribute='orthologs')
#' getProtein(id=c('YEAST00001','YEAST00002','YEAST00012'))
#' getProtein(id=c('YEAST00001','YEAST00002','YEAST00012'), attribute='gene_ontology')


getProtein <- function(id, attribute = NULL) {

    if (!is.null(attribute) && !(attribute %in% c("domains", "homeologs", "gene_ontology",
        "orthologs", "locus"))) {
        stop("You must provide a valid attribute.")
    }

    if (length(id) == 1) {
        url <- urlGenerator(endpoint = "protein", id = id, detail = attribute)
    } else if (length(id) > 1) {
        MAXLEN <- 100
        nr_req <- ceiling(length(id)/MAXLEN)
        req_nr <- rep(1:nr_req, each = MAXLEN)[1:length(id)]
        split_ids <- split(id, req_nr)
        body <- lapply(split_ids, FUN = function(x) {
            jsonlite::toJSON(list(ids = x, auto_unbox = TRUE))
        })

        url <- urlGenerator(endpoint = "protein", id = "bulk_retrieve")
        data <- requestFactory(url = url, body = body)
        # make all proteins objects
        data <- lapply(data, FUN = function(prot) {
            n <- names(prot)
            objectFactory(n, prot)
        })
        names(data) <- id

        if (is.null(attribute)) {
            return(data)
        } else {
            attribute_data <- lapply(data, function(x) {
                if (grepl("https://", x[[attribute]])) {
                  requestFactory(x[[attribute]])
                } else {
                  x[[attribute]]
                }
            })

            if (attribute == "locus") {
                # list of g range objects - needs to be merged into a single object
                attribute_data <- suppressWarnings(do.call("c", unlist(attribute_data,
                  use.names = FALSE)))
            } else {
                attribute_data <- plyr::rbind.fill(attribute_data[sapply(attribute_data,
                  function(x) class(x) == "data.frame")])
            }
            return(attribute_data)
        }
    }
    return(requestFactory(url))
}



