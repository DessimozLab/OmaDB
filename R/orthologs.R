#' An example orthologs object. 
#'
#' A dataframe containing information for the orthologs of protein YEAST00058.
#'
#' @format A dataframe object with 15 variables:
#' \describe{
#'   \item{entry_nr}{entry number of the ortholog}
#'   \item{omaid}{oma identifier of the ortholog}
#'   \item{canonicalid}{canonicalid of the ortholog}
#'   \item{sequence_md5}{sequence_md5 of the ortholog}
#'   \item{oma_group}{oma_group of the ortholog}
#'   \item{oma_hog_id}{hog id of the ortholog}
#'   \item{chromosome}{chromosomal location of the ortholog}
#'   \item{locus.start}{start locus of the ortholog}
#'   \item{locus.end}{end locus of the ortholog}
#'   \item{locus.strand}{locus strand of the ortholog}
#'   \item{is_main_isoform}{true/false}
#'   \item{rel_type}{relationship type of the ortholog to the gene}
#'   \item{distance}{ortholog distance}
#'   \item{score}{ortholog score}
#'   ...
#' }
#' @source \url{https://omabrowser.org/api/protein/YEAST00058/orthologs}
"orthologs"