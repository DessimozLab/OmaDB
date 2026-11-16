## Conversion of parsed JSON (nested lists as returned by
## jsonlite::parse_json(simplifyVector = FALSE)) into R objects.
##
## The rules are applied consistently, independent of the number of
## elements in a response:
##  * a JSON object becomes an "omadb_obj" (a named list),
##  * a JSON array of objects becomes a data.frame (nested objects are
##    flattened into columns named "<field>.<subfield>", arrays become list
##    columns),
##  * a JSON array of scalars becomes an atomic vector,
##  * protein and cDNA sequences become Biostrings objects and genomic loci
##    become GRanges objects.

.isRecord <- function(x) {
    is.list(x) && !is.null(names(x))
}

.isRecordList <- function(x) {
    is.list(x) && is.null(names(x)) && length(x) > 0L &&
        all(vapply(x, function(e) is.null(e) || .isRecord(e), logical(1)))
}

.isScalarList <- function(x) {
    is.list(x) && is.null(names(x)) &&
        all(vapply(x, .isScalarOrNull, logical(1)))
}

.isScalarOrNull <- function(v) {
    is.null(v) || (is.atomic(v) && length(v) == 1L)
}

.nullToNA <- function(vals) {
    vals[vapply(vals, is.null, logical(1))] <- list(NA)
    vals
}

## Convert any parsed JSON value into the corresponding R object.
.toR <- function(x) {
    if (.isRecord(x)) {
        return(.toObject(x))
    }
    if (is.list(x) && length(x) == 0L) {
        return(data.frame())
    }
    if (.isRecordList(x)) {
        return(.toDataFrame(x))
    }
    .convertValue(x)
}

.convertValue <- function(v) {
    if (is.null(v)) {
        return(NA)
    }
    if (is.atomic(v)) {
        return(v)
    }
    if (.isRecord(v)) {
        return(.toObject(v))
    }
    if (length(v) == 0L) {
        return(v)
    }
    if (.isScalarList(v)) {
        return(unlist(.nullToNA(v), use.names = FALSE))
    }
    if (.isRecordList(v)) {
        return(.toDataFrame(v))
    }
    lapply(v, .convertValue)
}

## Convert a JSON object into an omadb_obj.
.toObject <- function(rec) {
    fields <- names(rec)
    out <- lapply(fields, function(nm) .convertField(nm, rec[[nm]], rec))
    names(out) <- fields
    class(out) <- "omadb_obj"
    out
}

.convertField <- function(name, value, rec) {
    if (is.null(value)) {
        return(NA)
    }
    if (name == "sequence" && is.character(value)) {
        return(Biostrings::AAString(value))
    }
    if (name == "cdna" && is.character(value)) {
        value <- gsub("[^ACGTURYKMSWBDHVN.-]", "N", toupper(value))
        return(Biostrings::DNAString(value))
    }
    if (name == "locus" && is.list(value) && length(value) > 0L &&
            !is.null(rec[["omaid"]])) {
        exons <- if (.isRecord(value)) list(value) else value
        return(.lociToGRanges(rec[["omaid"]], rec[["chromosome"]],
            list(exons)))
    }
    .convertValue(value)
}

## Convert a list of JSON objects into a data.frame.
.toDataFrame <- function(records) {
    records <- lapply(records, function(r) if (is.null(r)) list() else r)
    flat <- lapply(records, .flattenRecord)
    cols <- unique(unlist(lapply(flat, names), use.names = FALSE))
    columns <- lapply(cols, function(col) {
        .simplifyColumn(lapply(flat, function(r) r[[col]]))
    })
    names(columns) <- cols
    .makeDataFrame(columns, length(records))
}

.makeDataFrame <- function(columns, nrow) {
    structure(columns, row.names = c(NA_integer_, -as.integer(nrow)),
        class = "data.frame")
}

.flattenRecord <- function(rec) {
    if (length(rec) == 0L) {
        return(list())
    }
    parts <- lapply(names(rec), function(nm) {
        v <- rec[[nm]]
        if (.isRecord(v)) {
            sub <- .flattenRecord(v)
            if (length(sub) > 0L) {
                names(sub) <- paste(nm, names(sub), sep = ".")
            }
            return(sub)
        }
        stats::setNames(list(v), nm)
    })
    do.call(c, parts)
}

.simplifyColumn <- function(vals) {
    if (all(vapply(vals, .isScalarOrNull, logical(1)))) {
        return(unlist(.nullToNA(vals), use.names = FALSE))
    }
    lapply(vals, .convertValue)
}

## Build a GRanges object from loci.
##
## `ids` and `chromosomes` are vectors with one element per protein, `loci`
## is a list with one element per protein, each being a list of exons
## (objects with start, end and strand).
.lociToGRanges <- function(ids, chromosomes, loci) {
    if (is.null(chromosomes)) {
        chromosomes <- rep(NA_character_, length(ids))
    }
    nExons <- lengths(loci)
    exons <- unlist(loci, recursive = FALSE)
    getField <- function(f) {
        unlist(.nullToNA(lapply(exons, function(e) e[[f]])),
            use.names = FALSE)
    }
    seqnames <- rep(as.character(chromosomes), nExons)
    seqnames[is.na(seqnames)] <- "unknown"
    GenomicRanges::GRanges(
        seqnames = seqnames,
        ranges = IRanges::IRanges(start = as.integer(getField("start")),
            end = as.integer(getField("end")),
            names = rep(as.character(ids), nExons)),
        strand = .toStrand(getField("strand")))
}

.toStrand <- function(x) {
    x <- as.character(x)
    out <- rep("*", length(x))
    out[x %in% c("1", "+")] <- "+"
    out[x %in% c("-1", "-")] <- "-"
    out
}
