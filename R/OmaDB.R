#' OmaDB: R wrapper for the OMA Browser REST API
#'
#' OmaDB is a wrapper for the REST API of the Orthologous MAtrix project
#' (OMA), a database for the inference of orthologs among complete genomes.
#' For more details on the OMA project, see \url{https://omabrowser.org/}.
#'
#' @section OmaDB functions:
#' The package contains a range of functions to query the database:
#' \itemize{
#'     \item \code{\link{getProtein}}, \code{\link{searchProtein}} and
#'     \code{\link{mapSequence}} to retrieve proteins, their orthologs,
#'     GO annotations, cross-references, \ldots
#'     \item \code{\link{getHOG}}, \code{\link{getHOGs}} and
#'     \code{\link{getOrthoXML}} for Hierarchical Orthologous Groups
#'     \item \code{\link{getOMAGroup}} for OMA Groups
#'     \item \code{\link{getGenome}}, \code{\link{getGenomePairs}} and
#'     \code{\link{getSharedAncestry}} for genomes and the orthologs
#'     between them
#'     \item \code{\link{getTaxonomy}} for the taxonomy
#'     \item \code{\link{getSynteny}} for (ancestral) synteny
#'     \item \code{\link{annotateSequence}} for a fast function prediction of
#'     a sequence
#' }
#'
#' The retrieved data is returned as data.frames or as objects of class
#' \code{\link{omadb_obj}}, using Bioconductor classes such as
#' \code{GRanges} and \code{AAString} where appropriate. Functions such as
#' \code{\link{getLocus}}, \code{\link{getTree}}, \code{\link{formatTopGO}}
#' and \code{\link{getTopGO}} convert the data into objects of the
#' GenomicRanges, ape and topGO packages.
#'
#' All responses from the OMA Browser are cached, see
#' \code{\link{setOmaCache}}.
#'
#' The functions are described in more detail in the package vignettes:
#' \code{browseVignettes("OmaDB")}.
#'
#' @keywords internal
#' @importFrom methods is new
#' @importFrom utils packageVersion URLencode head
#' @importFrom stats setNames
#' @importFrom Biostrings AAString DNAString
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
"_PACKAGE"
