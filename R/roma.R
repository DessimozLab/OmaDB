#' roma: A package for the orthology prediction data download from OMA database.
#'
#' roma is a wrapper for the REST API for the Orthologous MAtrix project (OMA) which is a  database for the inference of orthologs among complete genomes. For more details on the OMA project, see https://omabrowser.org/oma/home/.
#' 
#' 
#' @section  roma functions:
#' The package contains a range of functions that are used to query the database in an R friendly way. 
#'
#' 		- [formatTopGO](formatTopGO.Rd)
#'
#' 		- [getAnnotation](getAnnotation.Rd)
#'
#' 		- [getAttribute](getAttribute.Rd)
#'
#' 		- [getData](getData)
#'
#' 		- [getGenomeAlignment](getGenomeAlignment.Rd)
#'
#' 		- [getGRanges](getGRanges.Rd)
#'
#' 		- [getObjectAttributes](getObjectAttributes.Rd)
#'
#' 		- [getOntologies](getOntologies.Rd)
#'
#' 		- [getSequences](getSequences.Rd)
#'
#' 		- [getHOG](getHOG.Rd)
#'
#' 		- [getTaxonomy](getTaxonomy.Rd)
#'
#' 		- [getTopGO](getTopGO.Rd)
#'
#' 		- [getTree](getTree.Rd)
#'
#' 		- [getXref](getXref.Rd)
#'
#' 		- [mapSequence](mapSequence.Rd)
#'
#' 		- [resolveURL](resolveURL.Rd)
#'
#' The above functions are described in more detail in the package vignette's listed below:
#'
#'  [Get started with Roma] (library/roma/vignettes/roma.html)
#'
#'	[Exploring Hierarchical orthologous groups with roma] (library/roma/vignettes/exploring_hogs.html)
#'
#'	[Exploring Taxonomic trees with roma] (library/roma/vignettes/tree_visualisation.html)
#'
#'	[Sequence Analysis with roma] (library/roma/vignettes/sequence_mapping.html)
#'
"_PACKAGE"