parse <- function(json) jsonlite::parse_json(json, simplifyVector = FALSE)

test_that("a JSON object becomes an omadb_obj", {
    obj <- .toR(parse('{"a": 1, "b": "x", "c": null, "d": [1, 2]}'))
    expect_s3_class(obj, "omadb_obj")
    expect_identical(obj$a, 1L)
    expect_identical(obj$b, "x")
    expect_true(is.na(obj$c))
    expect_identical(obj$d, c(1L, 2L))
})

test_that("a list of objects is always a data.frame (issue #4)", {
    one <- .toR(parse('[{"id": "a", "GO_term": "GO:1", "score": 0.5}]'))
    expect_s3_class(one, "data.frame")
    expect_identical(nrow(one), 1L)
    expect_named(one, c("id", "GO_term", "score"))

    two <- .toR(parse(paste0('[{"id": "a", "GO_term": "GO:1"},',
        '{"id": "b", "GO_term": "GO:2"}]')))
    expect_s3_class(two, "data.frame")
    expect_identical(two$GO_term, c("GO:1", "GO:2"))

    none <- .toR(parse("[]"))
    expect_s3_class(none, "data.frame")
    expect_identical(nrow(none), 0L)
})

test_that("records are flattened and missing values become NA", {
    df <- .toR(parse(paste0(
        '[{"omaid": "A1", "species": {"code": "AAA", "taxon_id": 1},',
        '  "levels": ["x", "y"]},',
        ' {"omaid": "B1", "species": {"code": "BBB", "taxon_id": null},',
        '  "levels": [], "extra": true}]')))
    expect_named(df, c("omaid", "species.code", "species.taxon_id",
        "levels", "extra"))
    expect_identical(df$species.code, c("AAA", "BBB"))
    expect_identical(df$species.taxon_id, c(1L, NA))
    expect_identical(df$extra, c(NA, TRUE))
    expect_type(df$levels, "list")
    expect_identical(df$levels[[1]], c("x", "y"))
})

test_that("nested lists of records become data.frames", {
    obj <- .toR(parse(paste0('{"group_nr": 1, "members": [',
        '{"omaid": "A1"}, {"omaid": "B1"}]}')))
    expect_s3_class(obj$members, "data.frame")
    expect_identical(obj$members$omaid, c("A1", "B1"))
})

test_that("sequences and loci are converted to Bioconductor classes", {
    obj <- .toR(parse(paste0('{"omaid": "YEAST00001", "chromosome": "I",',
        '"locus": {"start": 10, "end": 20, "strand": -1},',
        '"sequence": "MKV", "cdna": "ATGAAXGTT"}')))
    expect_s4_class(obj$sequence, "AAString")
    expect_s4_class(obj$cdna, "DNAString")
    expect_identical(as.character(obj$cdna), "ATGAANGTT")
    expect_s4_class(obj$locus, "GRanges")
    expect_identical(as.character(GenomicRanges::seqnames(obj$locus)), "I")
    expect_identical(as.character(GenomicRanges::strand(obj$locus)), "-")
    expect_identical(names(obj$locus), "YEAST00001")

    exons <- .toR(parse(paste0('{"omaid": "X1", "chromosome": "2",',
        '"locus": [{"start": 1, "end": 5, "strand": "+"},',
        '{"start": 9, "end": 12, "strand": "+"}]}')))
    expect_length(exons$locus, 2L)
    expect_identical(GenomicRanges::start(exons$locus), c(1L, 9L))
})

test_that("as.data.frame converts an omadb_obj into a single row", {
    obj <- .toR(parse(paste0('{"omaid": "A1", "species": {"code": "AAA"},',
        '"levels": ["x", "y"]}')))
    df <- as.data.frame(obj)
    expect_s3_class(df, "data.frame")
    expect_identical(nrow(df), 1L)
    expect_named(df, c("omaid", "species.code", "levels"))
    expect_identical(df$levels[[1]], c("x", "y"))
})

test_that("print and getObjectAttributes describe objects", {
    obj <- .toR(parse(paste0('{"omaid": "A1", "orthologs": ',
        '"https://omabrowser.org/api/protein/A1/orthologs/"}')))
    expect_output(print(obj), "loaded on access")
    expect_identical(getObjectAttributes(obj),
        c(omaid = "character", orthologs = "url"))
})
