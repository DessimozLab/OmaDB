#' @title Deprecated functions in package \pkg{OmaDB}.
#' @description These functions are provided for compatibility with
#' older versions of \pkg{OmaDB} only, and will be defunct at the
#' next release.
#'
#' @details The following functions are deprecated and will be
#' made defunct; use the replacement indicated below:
#'
#' @name OmaDB-deprecated
#' @keywords internal
NULL



#' Get the CrossReferences in the OMA database for a pattern
#'
#' This function is should no longer be used. Use instead \code{\link{searchProtein}}.
#'
#' @param pattern the pattern to query the OMA database with - needs to be at least 3 characters long
#' @return a data.frame containing information on the cross references for a given pattern
#'
#' @templateVar fun getXref
#' @template template-depr_fun
NULL

#' @templateVar old getXref
#' @templateVar new searchProtein
#' @template template-depr_pkg
#'
#' @export
getXref <- function(pattern) {
    .Deprecated("searchProtein")
    return(searchProtein(pattern))
}

#' Get GO annotation for a sequence Function
#'
#' This function should no longer be used. Use instead \code{\link{annotateSequence}}.
#'
#' @param query the sequence to be annotated, it can be either a string or an AAString object from the Biostrings package
#' @return a data.frame containg the GO annotaition information linked to the query sequence
#'
#' @templateVar fun getAnnotation
#' @template template-depr_fun
NULL

#' @templateVar old getAnnotation
#' @templateVar new annotateSequence
#' @template template-depr_pkg
#' @export
getAnnotation <- function(query) {
    .Deprecated("annotateSequence")
    return(annotateSequence(query))
}

#' Get Whole Genome Alignment Function
#'
#' This function should no longer be used. Use instaed \code{\link{getGenomePairs}}.
#'
#' @param genome_id1 an identifier for the first genome, which can be either its taxon id or UniProt species code
#' @param genome_id2 an an identifier for the second genome, which can be either its taxon id or UniProt species code
#' @param chr1 the chromosome of interest for the first genome
#' @param chr2 the chromosome of interest for the second genome
#' @param rel_type the pairs relationship type
#' @param per_page the number of instances to be returned or 'all'. default is set to a 100.
#' @return a dataframe containing information about both the entries in the orthologous pair and their relationship
#'
#' @templateVar fun getGenomeAlignment
#' @template template-depr_fun
NULL

#' @templateVar old getGenomeAlignment
#' @templateVar new getGenomePairs
#' @template template-depr_pkg
#'
#' @export
getGenomeAlignment <- function(genome1, genome2, chr1 = NULL, chr2 = NULL, rel_type = NULL) {
    .Deprecated("getGenomePair")
    return(getGenomePairs(genome1, genome2, chr1, chr2, rel_type))
}


#' Get the Data Function
#'
#' The function to obtain the information available for a single entry in the datase.
#' This function should no longer be used. It has been divided into several functions:
#' Use the following functions instead.
#' \itemize{
#'   \item \code{\link{getProtein}} to obtain proteins (former \code{type='protein'})
#'   \item \code{\link{getGenome}} to obtain genomes (former \code{type='genome'})
#'   \item \code{\link{getOMAGroup}} to obtain genomes (former \code{type='group'})
#' }
#'
#' @param type the type for the entry to be returned - either protein, genome or group
#' @param id an identifier for the entry to be returned. For more information, see the 'Get started with OmaDB' vignette.
#' @param attribute an extra attribute
#' @return an object containing the JSON keys as attributes
#'
#' @templateVar fun getData
#' @template template-depr_fun
NULL

#' @templateVar old getData
#' @templateVar new getProtein}, \link{getGenome}, \link{getOMAGroup
#' @template template-depr_pkg
#'
#' @export
getData <- function(type, id, attribute = NULL) {
    if (type == "group") {
        .Deprecated("getOMAGroup")
        return(getOMAGroup(id, attribute))
    } else if (type == "genome") {
        .Deprecated("getGenome")
        return(getGenome(id, attribute))
    } else if (type == "protein") {
        .Deprecated("getProtein")
        return(getProtein(id, attribute))
    } else {
        .Deprecated("getProtein")
        return(NULL)
    }
}
