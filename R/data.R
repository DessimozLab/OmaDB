#' Example data retrieved from the OMA Browser
#'
#' Example responses of the OMA Browser REST API, as returned by the
#' functions of this package. They were retrieved from the OMA Browser
#' release All.May2026 with the script
#' \code{system.file("scripts", "make-data.R", package = "OmaDB")}.
#'
#' \describe{
#'     \item{protein}{\code{getProtein("YEAST00001")}: an
#'     \code{\link{omadb_obj}} with the data of a protein.}
#'     \item{orthologs}{\code{getProtein("YEAST00001", attribute =
#'     "orthologs")}: a data.frame with the orthologs of the protein.}
#'     \item{group}{\code{getOMAGroup("YEAST00058")}: an
#'     \code{\link{omadb_obj}} with the OMA Group of protein YEAST00058,
#'     including its members.}
#'     \item{hog}{\code{getHOG("YEAST00001")}: a data.frame with the HOG of
#'     protein YEAST00001.}
#'     \item{pairs}{\code{head(getGenomePairs("YEAST", "ASHGO", chr1 = "I"),
#'     100)}: a data.frame with orthologous pairs between the genomes of
#'     Saccharomyces cerevisiae and Ashbya gossypii.}
#'     \item{taxonomy}{\code{getTaxonomy(members = c("YEAST", "ASHGO",
#'     "SCHPO"))}: an \code{\link{omadb_obj}} with the taxonomy of three
#'     yeast species in newick format.}
#'     \item{xref}{\code{searchProtein("P53_RAT")}: a data.frame with the
#'     cross-references matching the query.}
#'     \item{sequence_map}{\code{mapSequence(seq)}: an \code{\link{omadb_obj}}
#'     with the proteins matching a query sequence.}
#'     \item{sequence_annotation}{\code{annotateSequence(seq)}: a data.frame
#'     with the predicted GO annotations of a query sequence.}
#' }
#'
#' @name OmaDB-data
#' @aliases protein orthologs group hog pairs taxonomy xref sequence_map
#'     sequence_annotation
#' @docType data
#' @format see the description of the corresponding function.
#' @source \url{https://omabrowser.org/api/}
#' @keywords datasets
#' @examples
#' data(protein, package = "OmaDB")
#' protein$omaid
#' data(pairs, package = "OmaDB")
#' table(pairs$rel_type)
NULL
