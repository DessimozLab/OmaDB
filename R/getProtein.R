.proteinAttributes <- c("orthologs", "hog_derived_orthologs", "homoeologs",
    "gene_ontology", "domains", "xref", "isoforms", "locus")

#' Retrieve a protein from the OMA Browser
#'
#' This function enables to retrieve information on one or several proteins
#' from the OMA Browser database.
#'
#' In its simplest form the function returns the base data of the query
#' protein. The query protein can be selected with any unique id, for example
#' with a UniProtKB accession (P12345), an OMA id (YEAST00012), or a RefSeq
#' id (NP_001226). To retrieve more than one protein, pass a vector of IDs.
#'
#' Non-scalar properties of proteins such as their domains, GO annotations,
#' orthologs or homoeologs are loaded upon accessing them, or - if you only
#' need this information - you can set the attribute parameter to the
#' property name and retrieve this information directly:
#' \describe{
#'     \item{orthologs}{the orthologs of the protein, optionally restricted to
#'     a type of relationship with \code{rel_type}}
#'     \item{hog_derived_orthologs}{the orthologs of the protein derived from
#'     its hierarchical orthologous groups}
#'     \item{homoeologs}{the homoeologs of the protein (only for proteins of
#'     polyploid genomes)}
#'     \item{gene_ontology}{the Gene Ontology annotations of the protein}
#'     \item{domains}{the domain annotations of the protein}
#'     \item{xref}{the cross-references of the protein to other databases}
#'     \item{isoforms}{the splicing isoforms of the protein's gene}
#'     \item{locus}{the genomic locus of the protein as a GRanges object}
#' }
#'
#' @seealso For non-unique IDs or partial ID lookup, use
#'     \code{\link{searchProtein}} instead.
#'
#' @param id Identifier(s) for the entry or entries to be returned. A
#'     character string for a single entry or a vector for several entries.
#' @param attribute Instead of the protein, return the attribute property of
#'     the protein. See details for the available attributes. 'homeologs' is
#'     accepted as an alias for 'homoeologs'.
#' @param rel_type For \code{attribute = "orthologs"}, restrict the orthologs
#'     to a type of relationship, i.e. '1:1', '1:n', 'm:1' or 'm:n'.
#' @return For a single id, an \code{\link{omadb_obj}} with the protein's
#'     attributes. For several ids, a named list of such objects (with
#'     \code{NULL} for ids that are not found). If \code{attribute} is given,
#'     a data.frame with the requested data (combined for all ids, with a
#'     column \code{query_id}), except for 'locus' which returns a GRanges
#'     object and 'domains' which returns an \code{omadb_obj} for a single id
#'     and a list of them for several ids.
#' @export
#' @examples
#' getProtein(id = "YEAST00001")
#' head(getProtein(id = "YEAST00001", attribute = "orthologs"))
#' prots <- getProtein(id = c("YEAST00001", "YEAST00002", "YEAST00012"))
#' names(prots)
#' go <- getProtein(id = c("YEAST00001", "YEAST00002"),
#'     attribute = "gene_ontology")
#' table(go$query_id)
getProtein <- function(id, attribute = NULL, rel_type = NULL) {
    if (missing(id) || length(id) == 0L) {
        stop("You must provide at least one protein id.")
    }
    if (!is.null(attribute)) {
        if (identical(attribute, "homeologs")) {
            attribute <- "homoeologs"
        }
        attribute <- match.arg(attribute, .proteinAttributes)
    }
    .checkRelType(rel_type, c("1:1", "1:n", "m:1", "m:n"))
    if (length(id) == 1L) {
        return(.proteinData(id, attribute, rel_type))
    }
    if (is.null(attribute)) {
        return(.bulkRetrieveProteins(id))
    }
    if (attribute == "domains") {
        res <- lapply(id, .proteinData, attribute = attribute)
        names(res) <- id
        return(res)
    }
    parts <- lapply(id, .proteinData, attribute = attribute,
        rel_type = rel_type)
    if (attribute == "locus") {
        return(do.call(c, unname(parts)))
    }
    .rbindWithQueryId(parts, id)
}

.proteinData <- function(id, attribute = NULL, rel_type = NULL) {
    .checkScalarId(id)
    if (identical(attribute, "locus")) {
        prot <- .toR(.omaGet(.apiUrl("protein", id)))
        return(.subset2(prot, "locus"))
    }
    query <- if (identical(attribute, "orthologs")) {
        list(rel_type = rel_type)
    } else {
        list()
    }
    url <- .apiUrl("protein", id, attribute, query = query)
    .toR(.omaGet(url))
}

## Retrieve many proteins at once using the bulk_retrieve endpoint, which
## accepts at most 100 ids per request.
.bulkRetrieveProteins <- function(ids) {
    chunkSize <- 100L
    chunks <- split(as.character(ids),
        ceiling(seq_along(ids) / chunkSize))
    url <- .apiUrl("protein", "bulk_retrieve")
    results <- lapply(chunks, function(chunk) {
        body <- jsonlite::toJSON(list(ids = chunk))
        .omaPost(url, body = as.character(body))
    })
    results <- unlist(unname(results), recursive = FALSE)
    out <- lapply(results, function(r) {
        if (is.null(r$target)) NULL else .toObject(r$target)
    })
    names(out) <- vapply(results, function(r) {
        as.character(r$query_id)
    }, character(1))
    out[as.character(ids)]
}

.checkScalarId <- function(id) {
    if (length(id) != 1L || is.na(id) || !nzchar(as.character(id))) {
        stop("You must provide a single, valid id.")
    }
    invisible(TRUE)
}

.checkRelType <- function(rel_type, valid) {
    if (!is.null(rel_type) && !(length(rel_type) == 1L &&
            rel_type %in% valid)) {
        stop("rel_type must be one of ",
            paste(sQuote(valid, FALSE), collapse = ", "), ".")
    }
    invisible(TRUE)
}

## Combine data.frames retrieved for several query ids, adding a column
## with the query id. Results that are empty are skipped.
.rbindWithQueryId <- function(parts, ids) {
    parts <- lapply(seq_along(parts), function(i) {
        df <- parts[[i]]
        if (!is.data.frame(df) || nrow(df) == 0L) {
            return(NULL)
        }
        cbind(query_id = rep(as.character(ids[[i]]), nrow(df)), df,
            stringsAsFactors = FALSE)
    })
    parts <- Filter(Negate(is.null), parts)
    if (length(parts) == 0L) {
        return(data.frame())
    }
    .rbindFill(parts)
}

## rbind data.frames with possibly differing columns, filling missing
## columns with NA.
.rbindFill <- function(dfs) {
    cols <- unique(unlist(lapply(dfs, names), use.names = FALSE))
    dfs <- lapply(dfs, function(df) {
        for (col in setdiff(cols, names(df))) {
            df[[col]] <- rep(NA, nrow(df))
        }
        df[cols]
    })
    out <- do.call(rbind, dfs)
    rownames(out) <- NULL
    out
}
