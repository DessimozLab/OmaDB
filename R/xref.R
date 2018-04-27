#' An example xref object. 
#'
#'
#' @format A dataframe with 8 variables:
#' \describe{
#'   \item{xref}{cross reference}
#'   \item{source}{source of the cross reference}
#'   \item{entry_nr}{oma database entry number}
#'   \item{oma_id}{oma id of the cross reference}
#'   \item{genome.code}{genome_id of the cross reference}
#'   \item{genome.taxon_id}{taxon_id of the cross reference}
#'   \item{genome.species}{species of the cross reference}
#'   \item{genome.genome_url}{genome url pointer of the cross reference}
#'   ...
#' }
#' @source \url{https://omabrowser.org/api/xref/?search=MAL}
"xref"