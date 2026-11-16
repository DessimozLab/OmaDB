#' Map a protein sequence to the OMA Browser
#'
#' The function identifies the protein(s) in the OMA Browser that match a
#' query sequence, either exactly or approximately.
#'
#' @param query the sequence to be searched, it can be either a string or an
#'     AAString object from the Biostrings package
#' @param search argument to choose search strategy. Can be set to 'exact',
#'     'approximate' or 'mixed'. Defaults to 'mixed', meaning first tries to
#'     find exact match. If no target can be found, uses approximate search
#'     strategy to identify query sequence in database.
#' @param full_length a logical indicating whether or not for exact matches,
#'     the query sequence must be matching the full target sequence. By
#'     default, a partial exact match is also reported as exact match.
#' @return an \code{\link{omadb_obj}} with the attributes \code{query},
#'     \code{identified_by} (the type of match) and \code{targets}, a
#'     data.frame with the matching proteins.
#' @export
#' @examples
#' seq <- paste0("MNDPSLLGYPNVGPQQQQQQQQQQHAGLLGKGTPNALQQQLHMNQLTGIPPPGLMN",
#'     "NSDVHTSSNNNSRQLLDQLANGNANMLNMNMDNNNNNNNNNNNNNNNGGGSGVMMNASTAAVNSIGM")
#' res <- mapSequence(query = seq)
#' res$identified_by
#' res$targets[, c("omaid", "canonicalid")]
mapSequence <- function(query, search = NULL, full_length = FALSE) {
    if (missing(query)) {
        stop("You must provide a sequence to query.")
    }
    if (!is.null(search)) {
        search <- match.arg(tolower(search),
            c("mixed", "exact", "approximate"))
    }
    query <- .sequenceArg(query)
    url <- .apiUrl("sequence", query = list(query = query,
        full_length = full_length, search = search))
    .toR(.omaGet(url))
}

.sequenceArg <- function(query) {
    if (is(query, "XString")) {
        query <- as.character(query)
    }
    if (!is.character(query) || length(query) != 1L || !nzchar(query)) {
        stop("query must be a single sequence.")
    }
    query
}
