#' Retrieve a genome from the OMA Browser database
#'
#' This function obtains the basic information for one specific genome
#' available on the OMA Browser, or - if no id is provided - a data.frame
#' with all available genomes.
#'
#' Ids can be either the scientific name of a species, the NCBI taxonomy id
#' or the UniProtKB mnemonic species code.
#'
#' The optional argument attribute can be used to directly load the proteins
#' (all proteins, including alternative splicing isoforms) or the genes
#' (only the main isoform of each gene) of the genome. Alternatively, you can
#' access the \code{proteins} attribute of the result which will
#' transparently load the proteins from the OMA Browser.
#'
#' @param id A genome identifier. By default, all available genomes will be
#'     returned.
#' @param attribute An extra attribute to be returned ('proteins' or
#'     'genes').
#' @return If no id is given, a data.frame with all genomes. Otherwise an
#'     \code{\link{omadb_obj}} with the genome's attributes, or - if
#'     \code{attribute} is given - a data.frame with the proteins or genes of
#'     the genome.
#' @export
#' @examples
#' genomes <- getGenome()
#' head(genomes[, c("code", "taxon_id", "species")])
#' getGenome(id = "HUMAN")
#' getGenome(id = 9606)
#' prots <- getGenome(id = "ASHGO", attribute = "proteins")
#' nrow(prots)
getGenome <- function(id = NULL, attribute = NULL) {
    if (!is.null(attribute)) {
        attribute <- match.arg(attribute, c("proteins", "genes"))
        if (is.null(id)) {
            stop("You must provide a genome id to retrieve its ",
                attribute, ".")
        }
    }
    if (is.null(id)) {
        return(.toR(.omaGetList(.apiUrl("genome"))))
    }
    .checkScalarId(id)
    if (is.null(attribute)) {
        return(.toR(.omaGet(.apiUrl("genome", id))))
    }
    .toR(.omaGetList(.apiUrl("genome", id, attribute)))
}
