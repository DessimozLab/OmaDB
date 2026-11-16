#' Objects returned by the OmaDB package
#'
#' Functions that retrieve a single entity from the OMA Browser (e.g. a
#' protein, a genome, an OMA group or a HOG) return an object of class
#' \code{omadb_obj}. It is a named list with one element per attribute of the
#' entity. Attributes that hold a list of entities (e.g. the members of an
#' OMA group) are data.frames, protein and cDNA sequences are
#' \code{\link[Biostrings]{AAString}} and \code{\link[Biostrings]{DNAString}}
#' objects and genomic loci are \code{\link[GenomicRanges]{GRanges}} objects.
#'
#' Some attributes are not included in the response from the OMA Browser,
#' but are given as a url pointing to the API endpoint that provides the
#' data. Accessing such an attribute with \code{$} (or
#' \code{\link{getAttribute}}) loads the data transparently. The responses
#' are cached (see \code{\link{setOmaCache}}), so accessing an attribute
#' repeatedly does not query the server again. Use \code{[[} to access the
#' url itself.
#'
#' \code{as.data.frame} converts an object into a data.frame with a single
#' row. Nested attributes are flattened into columns named
#' \code{<attribute>.<subattribute>}; attributes that are not scalars
#' (vectors, data.frames, sequences, \ldots) become list columns.
#'
#' @param x an object of class \code{omadb_obj}
#' @param name the name of the attribute
#' @param row.names,optional,... passed on to/ignored as for
#'     \code{\link[base]{as.data.frame}}
#' @return \code{$} returns the value of the attribute, loading it from the
#'     OMA Browser if necessary. \code{as.data.frame} returns a data.frame
#'     with one row. \code{print} returns \code{x} invisibly.
#' @name omadb_obj
#' @aliases omadb_obj-class
#' @examples
#' prot <- getProtein("YEAST00001")
#' prot$omaid
#' ## the orthologs are loaded on access
#' head(prot$orthologs[, c("omaid", "rel_type")])
#' as.data.frame(prot)[, c("omaid", "species.code", "oma_group")]
NULL

#' @rdname omadb_obj
#' @export
"$.omadb_obj" <- function(x, name) {
    value <- .subset2(x, name)
    if (.isApiUrl(value)) {
        value <- resolveURL(value)
    }
    value
}

#' @rdname omadb_obj
#' @export
as.data.frame.omadb_obj <- function(x, row.names = NULL, optional = FALSE,
        ...) {
    flat <- .flattenObject(x)
    columns <- lapply(flat, function(v) {
        if (is.atomic(v) && length(v) == 1L) v else list(v)
    })
    df <- .makeDataFrame(columns, 1L)
    if (!is.null(row.names)) {
        row.names(df) <- row.names
    }
    df
}

.flattenObject <- function(x) {
    parts <- lapply(names(x), function(nm) {
        v <- .subset2(x, nm)
        if (is(v, "omadb_obj")) {
            sub <- .flattenObject(v)
            if (length(sub) > 0L) {
                names(sub) <- paste(nm, names(sub), sep = ".")
            }
            return(sub)
        }
        stats::setNames(list(v), nm)
    })
    do.call(c, parts)
}

#' @rdname omadb_obj
#' @exportS3Method base::print
print.omadb_obj <- function(x, ...) {
    cat("OMA Browser object with", length(x), "attributes:\n")
    width <- max(nchar(names(x)), 0L)
    for (nm in names(x)) {
        cat(" $", formatC(nm, width = -width), ":",
            .describeValue(.subset2(x, nm)), "\n")
    }
    invisible(x)
}

.describeValue <- function(v) {
    if (.isApiUrl(v)) {
        return("<loaded on access>")
    }
    if (is.data.frame(v)) {
        return(sprintf("data.frame [%d x %d]", nrow(v), ncol(v)))
    }
    if (is(v, "omadb_obj")) {
        return(paste0("object (", paste(names(v), collapse = ", "), ")"))
    }
    if (is(v, "XString")) {
        return(sprintf("%s of length %d", class(v)[1], length(v)))
    }
    if (is(v, "GRanges")) {
        return(sprintf("GRanges with %d range(s)", length(v)))
    }
    if (is.atomic(v)) {
        txt <- paste(utils::head(format(v), 3L), collapse = " ")
        if (length(v) > 3L) {
            txt <- paste(txt, "...")
        }
        if (nchar(txt) > 60L) {
            txt <- paste0(substr(txt, 1L, 57L), "...")
        }
        return(sprintf("%s %s", class(v)[1], txt))
    }
    sprintf("%s of length %d", class(v)[1], length(v))
}

#' Get the attributes of an object
#'
#' Lists the attributes of an object returned by the package and their
#' data types. Attributes of type "url" are loaded from the OMA Browser when
#' they are accessed.
#'
#' @param obj an object of class \code{\link{omadb_obj}}
#' @return a named character vector with the data type of each attribute
#' @export
#' @examples
#' getObjectAttributes(getOMAGroup(id = "YEAST58"))
getObjectAttributes <- function(obj) {
    vapply(names(obj), function(nm) {
        v <- .subset2(obj, nm)
        if (.isApiUrl(v)) "url" else class(v)[1]
    }, character(1))
}

#' Get the value of an attribute of an object
#'
#' Returns the value of an attribute of an object returned by the package.
#' If the attribute is a url pointing to the OMA Browser API, the data
#' behind the url is loaded and returned. This is equivalent to
#' \code{obj$attribute}.
#'
#' @param obj an object of class \code{\link{omadb_obj}}
#' @param attribute the name of the attribute
#' @return the value of the attribute
#' @export
#' @examples
#' members <- getAttribute(getOMAGroup(id = "YEAST58"), "members")
#' head(members$omaid)
getAttribute <- function(obj, attribute) {
    stopifnot(is.character(attribute), length(attribute) == 1L)
    value <- .subset2(obj, attribute)
    if (.isApiUrl(value)) {
        value <- resolveURL(value)
    }
    value
}

#' Load data for a given url from the OMA Browser API
#'
#' Loads the data behind a url of the OMA Browser API. This is usually not
#' needed, as attributes containing a url are loaded automatically when
#' accessed with \code{$}. It is useful for urls that are stored in a
#' data.frame, e.g. the \code{members_url} column of the result of
#' \code{\link{getHOGs}}.
#'
#' The type of the result depends only on the type of data provided by the
#' endpoint: endpoints that provide a single entity return an
#' \code{\link{omadb_obj}}, endpoints that provide a list of entities always
#' return a data.frame (also if the list contains a single or no element).
#' Paginated endpoints are loaded completely.
#'
#' @param url the url of interest
#' @return an \code{\link{omadb_obj}} or a data.frame
#' @export
#' @examples
#' go <- resolveURL(paste0("https://omabrowser.org/api/protein/YEAST00001/",
#'     "gene_ontology/"))
#' head(go)
resolveURL <- function(url) {
    stopifnot(is.character(url), length(url) == 1L, !is.na(url))
    .toR(.omaGet(url))
}
