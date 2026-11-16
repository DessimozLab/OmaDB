#' Get the topGO Object
#'
#' The function creates a \code{topGOdata} object containing the GO
#' annotations for the given genes, which can be used to perform a GO
#' enrichment analysis with the topGO package.
#'
#' The topGO package must be attached with \code{library(topGO)} before
#' calling this function, as topGO requires it to find its GO term
#' environments.
#'
#' @param annotations list of GO annotations obtained from
#'     \code{\link{formatTopGO}}
#' @param format Format of the annotations - either 'GO2geneID' or
#'     'geneID2GO'
#' @param foregroundGenes identifiers of the foreground genes, i.e. the genes
#'     of interest. Either a character vector or a data.frame with a column
#'     \code{omaid}.
#' @param ontology The ontology for which the enrichment should be done
#'     ('BP', 'MF' or 'CC'). This parameter is passed directly to the
#'     topGOdata constructor.
#' @return a \code{topGOdata} object
#' @export
#' @import topGO
#' @importFrom methods new is
#' @examples
#' library(topGO)
#' members <- getHOG("HOG:F0411479", level = "Saccharomycetaceae",
#'     members = TRUE)
#' go <- getProtein(members$omaid, attribute = "gene_ontology")
#' annotations <- formatTopGO(go, format = "geneID2GO")
#' GOdata <- getTopGO(annotations, foregroundGenes = head(members$omaid, 3),
#'     format = "geneID2GO", ontology = "BP")
#' GOdata
getTopGO <- function(annotations, format, foregroundGenes, ontology) {
    if (missing(ontology)) {
        stop("You must provide the ontology for which the enrichment ",
            "should be done.")
    }
    if (missing(annotations) || !is.list(annotations)) {
        stop("You must provide a valid list of annotations.")
    }
    if (missing(format) || !(format %in% c("geneID2GO", "GO2geneID"))) {
        stop("You must provide a valid annotations format.")
    }
    if (missing(foregroundGenes)) {
        stop("You must provide a valid list of genes of interest.")
    }
    if (!("package:topGO" %in% search())) {
        stop("The topGO package must be attached: call library(topGO) ",
            "before using getTopGO().")
    }
    if (is.data.frame(foregroundGenes)) {
        foregroundGenes <- foregroundGenes[["omaid"]]
    }
    foregroundGenes <- as.character(unlist(foregroundGenes))
    geneNames <- if (format == "geneID2GO") {
        names(annotations)
    } else {
        unique(unlist(annotations, use.names = FALSE))
    }
    geneList <- factor(as.integer(geneNames %in% foregroundGenes),
        levels = c(0L, 1L))
    names(geneList) <- geneNames
    if (format == "geneID2GO") {
        new("topGOdata", ontology = ontology, allGenes = geneList,
            annot = topGO::annFUN.gene2GO, gene2GO = annotations)
    } else {
        new("topGOdata", ontology = ontology, allGenes = geneList,
            annot = topGO::annFUN.GO2genes, GO2genes = annotations)
    }
}
