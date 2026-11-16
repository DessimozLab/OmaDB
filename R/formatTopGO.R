#' Format the GO annotations data
#'
#' The function creates a list of GO annotations that is compatible with
#' topGO from proteins or GO annotations retrieved from the OMA Browser.
#'
#' @param geneList either a list of protein objects (as returned by
#'     \code{\link{getProtein}} for several ids), or a data.frame of GO
#'     annotations, e.g. as returned by \code{getProtein(ids, attribute =
#'     "gene_ontology")} or \code{getHOG(id, attribute = "gene_ontology")}.
#'     For a data.frame, the genes are identified by the column
#'     \code{query_id} if present, \code{id} otherwise.
#' @param format format for the data to be returned in - either 'GO2geneID'
#'     or 'geneID2GO'
#' @return a named list containing the GO2geneID or geneID2GO information
#' @export
#' @examples
#' geneList <- getProtein(id = c("YEAST00001", "YEAST00003"))
#' annotations <- formatTopGO(geneList, format = "geneID2GO")
#' lengths(annotations)
formatTopGO <- function(geneList, format) {
    if (missing(format) || !(format %in% c("GO2geneID", "geneID2GO"))) {
        stop("Invalid format. Must be either 'GO2geneID' or 'geneID2GO'")
    }
    if (is.data.frame(geneList)) {
        geneID2GO <- .goFromDataFrame(geneList)
    } else if (is.list(geneList)) {
        geneList <- Filter(Negate(is.null), geneList)
        geneID2GO <- lapply(geneList, function(protein) {
            go <- protein$gene_ontology
            if (is.data.frame(go) && "GO_term" %in% names(go)) {
                unique(go$GO_term)
            } else {
                character(0)
            }
        })
        names(geneID2GO) <- vapply(geneList, function(protein) {
            as.character(.subset2(protein, "omaid"))
        }, character(1))
    } else {
        stop("geneList must be a list of proteins or a data.frame.")
    }
    if (format == "geneID2GO") {
        return(geneID2GO)
    }
    topGO::inverseList(geneID2GO)
}

.goFromDataFrame <- function(df) {
    idCol <- intersect(c("query_id", "id"), names(df))
    if (length(idCol) == 0L || !("GO_term" %in% names(df))) {
        stop("The data.frame must have the columns GO_term and ",
            "query_id or id.")
    }
    genes <- as.character(df[[idCol[1]]])
    lapply(split(as.character(df$GO_term), factor(genes, unique(genes))),
        unique)
}
