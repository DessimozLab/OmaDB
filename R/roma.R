#' roma: A package for the orthology prediction data download from OMA database.
#'
#' roma is a wrapper for the REST API for the Orthologous MAtrix project (OMA) which is a  database for the inference of orthologs among complete genomes. For more details on the OMA project, see \url{https://omabrowser.org/oma/home/}.
#' 
#' 
#' @section  roma functions:
#' The package contains a range of functions that are used to query the database. Some of the main functions are listed below:
#'
#' \itemize{
#'   \item getData()
#'   \item getHOG()
#'   \item getGenomeAlignment()
#'   \item getTaxonomy()
#'   \item mapSequence()
#'   \item getAnnotation()
#'   \item getXref()
#' } 
#' 	
#' In addition to these, roma features a range of functions that are used to format the retrieved data into some commonly used Bioconductor objects using packages such as GenomicRanges, Biostrings, topGO and ggtree.	Some of them are listed below:
#'
#' \itemize{
#'   \item formatTopGO()
#'   \item getGRanges()
#'   \item getTree()
#' } 
#'
#' The above functions are described in more detail in the package vignette's listed below:
#'
#'
#' \itemize{
#'   \item Get started with Roma
#'   \item Exploring Hierarchical orthologous groups with roma}
#'   \item Exploring Taxonomic trees with roma
#'   \item Sequence Analysis with roma
#' }  
#'

"_PACKAGE"