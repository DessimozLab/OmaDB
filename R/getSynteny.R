#' Retrieve synteny information from the OMA Browser
#'
#' The function retrieves the (ancestral) synteny graph, i.e. the gene
#' order, at a given taxonomic level. Nodes of the graph are genes (HOGs for
#' ancestral levels, proteins for extant genomes), edges connect genes that
#' are (inferred to be) neighbours on a chromosome.
#'
#' If \code{id} is given, only the neighbourhood of the given HOG or
#' protein is returned, otherwise the synteny graph of the whole (ancestral)
#' genome.
#'
#' The result can easily be converted into a graph object, e.g. with
#' \code{igraph::graph_from_data_frame(res$links, directed = FALSE)}. Note
#' that for a neighbourhood, the links of the outermost nodes may point to
#' nodes that are not included in \code{nodes}.
#'
#' @param level the taxonomic level of interest. Ancestral levels accept the
#'     scientific name or the NCBI taxon id, extant genomes also accept the
#'     UniProt species code.
#' @param id optional, a HOG ID (for ancestral levels) or protein id (for
#'     extant genomes) to retrieve its neighbourhood.
#' @param evidence the minimum evidence of the edges of the ancestral
#'     synteny graph: 'linearized' (default), 'parsimonious' or 'any'.
#' @param context the size of the neighbourhood around \code{id}, in number
#'     of edges.
#' @param break_circular_contigs logical, whether to break circular
#'     ancestral contigs at their weakest edge. Only used if evidence is
#'     'linearized'.
#' @return a list with two data.frames: \code{nodes} with one row per gene
#'     (the first column \code{id} identifies the gene) and \code{links} with
#'     one row per edge (the first columns \code{source} and \code{target},
#'     followed by \code{weight} and further information on the edge).
#' @export
#' @examples
#' syn <- getSynteny(level = "Saccharomycetaceae", id = "HOG:F0411479",
#'     context = 2)
#' syn$nodes
#' syn$links[, c("source", "target", "weight")]
getSynteny <- function(level, id = NULL,
        evidence = c("linearized", "parsimonious", "any"), context = 2L,
        break_circular_contigs = TRUE) {
    if (missing(level)) {
        stop("You must provide a taxonomic level.")
    }
    .checkScalarId(level)
    evidence <- match.arg(evidence)
    query <- list(level = level, evidence = evidence,
        break_circular_contigs = break_circular_contigs)
    if (!is.null(id)) {
        .checkScalarId(id)
        query$context <- as.integer(context)
    }
    res <- .omaGet(.apiUrl("synteny", id, query = query))
    if (is.null(id)) {
        ## the whole genome is returned as a list of contigs
        contigs <- lapply(res, .syntenyGraph)
        return(list(
            nodes = .rbindFill(lapply(contigs, `[[`, "nodes")),
            links = .rbindFill(lapply(contigs, `[[`, "links"))))
    }
    .syntenyGraph(res)
}

.syntenyGraph <- function(graph) {
    nodes <- .toDataFrame(graph$nodes)
    if ("id" %in% names(nodes)) {
        nodes <- nodes[c("id", setdiff(names(nodes), "id"))]
    }
    links <- .toDataFrame(graph$links)
    first <- intersect(c("source", "target"), names(links))
    links <- links[c(first, setdiff(names(links), first))]
    list(nodes = nodes, links = links)
}
