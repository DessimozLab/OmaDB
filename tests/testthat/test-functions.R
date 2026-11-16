## Tests of exported functions against a fake API, for cases that are not
## covered by the recorded responses.

test_that("getHOGs returns all HOGs at a level", {
    calls <- local_fake_api(stats::setNames(list(
        paste0('[{"hog_id": "HOG:1", "level": "Fungi", ',
            '"completeness_score": 1},',
            '{"hog_id": "HOG:2", "level": "Fungi", ',
            '"completeness_score": 0.5}]'),
        '[{"hog_id": "HOG:1", "event": "duplicated"}]'),
        c(api("hog/?level=Fungi&per_page=5000"),
            api("hog/?level=Fungi&compare_with=Eukaryota&per_page=5000"))))
    hogs <- getHOGs("Fungi")
    expect_s3_class(hogs, "data.frame")
    expect_identical(hogs$hog_id, c("HOG:1", "HOG:2"))
    events <- getHOGs("Fungi", compare_with = "Eukaryota")
    expect_identical(events$event, "duplicated")
})

test_that("getSynteny combines the contigs of a whole genome", {
    url <- api(paste0("synteny/?level=Fungi&evidence=any&",
        "break_circular_contigs=false"))
    local_fake_api(stats::setNames(list(paste0(
        '[{"nodes": [{"nr_members": 3, "id": "A"}, {"id": "B"}],',
        '  "links": [{"weight": 2, "source": "A", "target": "B"}]},',
        ' {"nodes": [{"id": "C"}], "links": []}]')), url))
    syn <- getSynteny("Fungi", evidence = "any",
        break_circular_contigs = FALSE)
    expect_identical(syn$nodes$id, c("A", "B", "C"))
    expect_identical(syn$nodes$nr_members, c(3L, NA, NA))
    expect_identical(names(syn$links), c("source", "target", "weight"))
})

test_that("getLocus accepts data.frames, lists of them and objects", {
    df <- data.frame(omaid = c("A1", "B1"), chromosome = c("1", NA),
        locus.start = c(1L, 10L), locus.end = c(5L, 20L),
        locus.strand = c(1L, -1L))
    gr <- getLocus(df)
    expect_s4_class(gr, "GRanges")
    expect_identical(as.character(GenomicRanges::seqnames(gr)),
        c("1", "unknown"))
    expect_identical(as.character(GenomicRanges::strand(gr)), c("+", "-"))
    expect_identical(getLocus(list(df[1, ], df[2, ])), gr)
    expect_error(getLocus(data.frame(omaid = "A")), "must contain")
    expect_error(getLocus(list(1, 2)), "must be ids")

    obj <- .toR(jsonlite::parse_json(paste0('{"omaid": "A1",',
        '"locus": {"start": 1, "end": 5, "strand": 1}}')))
    expect_identical(getLocus(obj), obj[["locus"]])
})

test_that("getTree reads newick strings with quoted labels", {
    tree <- getTree('("A b (strain 1)",C)"Root level";')
    expect_identical(tree$tip.label, c("A b (strain 1)", "C"))
    expect_identical(tree$node.label, "Root level")
    expect_identical(getTree("(A,B)C;")$tip.label, c("A", "B"))
    expect_error(getTree(1), "newick string")
})

test_that("getTopGO creates a topGOdata object", {
    annotations <- list(G1 = c("GO:0006281", "GO:0006974"),
        G2 = "GO:0006281", G3 = "GO:0008150")
    expect_error(getTopGO(annotations, format = "geneID2GO",
        foregroundGenes = "G1"), "ontology")
    expect_error(getTopGO("x", format = "geneID2GO", foregroundGenes = "G1",
        ontology = "BP"), "list of annotations")
    expect_error(getTopGO(annotations, format = "foo",
        foregroundGenes = "G1", ontology = "BP"), "annotations format")
    expect_error(getTopGO(annotations, format = "geneID2GO",
        ontology = "BP"), "genes of interest")

    suppressMessages(withr::local_package("topGO"))
    GOdata <- suppressMessages(getTopGO(annotations, format = "geneID2GO",
        foregroundGenes = data.frame(omaid = "G1"), ontology = "BP"))
    expect_s4_class(GOdata, "topGOdata")
    expect_identical(sort(topGO::sigGenes(GOdata)), "G1")

    go2genes <- topGO::inverseList(annotations)
    GOdata <- suppressMessages(getTopGO(go2genes, format = "GO2geneID",
        foregroundGenes = "G1", ontology = "BP"))
    expect_s4_class(GOdata, "topGOdata")
})

test_that("formatTopGO validates its input", {
    expect_error(formatTopGO(list(), format = "foo"), "Invalid format")
    expect_error(formatTopGO("x", format = "geneID2GO"), "list of proteins")
    expect_error(formatTopGO(data.frame(a = 1), format = "geneID2GO"),
        "GO_term")
    hogGO <- data.frame(id = c("H1", "H1", "H2"),
        GO_term = c("GO:1", "GO:2", "GO:1"))
    expect_identical(formatTopGO(hogGO, format = "geneID2GO"),
        list(H1 = c("GO:1", "GO:2"), H2 = "GO:1"))
})

test_that("print describes all kinds of attributes", {
    obj <- .toR(jsonlite::parse_json(paste0(
        '{"omaid": "A1", "n": [1, 2, 3, 4],',
        '"members": [{"a": 1}], "species": {"code": "AAA"}, "nested": [[1]],',
        '"sequence": "MKV", "chromosome": "1",',
        '"locus": {"start": 1, "end": 5, "strand": 1},',
        '"description": "', strrep("x", 80), '"}')))
    out <- capture.output(print(obj))
    expect_match(out[1], "9 attributes")
    expect_true(any(grepl("data.frame \\[1 x 1\\]", out)))
    expect_true(any(grepl("object \\(code\\)", out)))
    expect_true(any(grepl("AAString of length 3", out)))
    expect_true(any(grepl("GRanges with 1 range", out)))
    expect_true(any(grepl("integer 1 2 3 \\.\\.\\.", out)))
    expect_true(any(grepl("list of length 1", out)))
    expect_true(any(grepl("xxx\\.\\.\\.", out)))
})

test_that("functions validate ids", {
    calls <- local_fake_api(list())
    expect_error(getOMAGroup(), "valid group id")
    expect_error(getOMAGroup(c("A", "B")), "single, valid id")
    expect_error(getSharedAncestry("A"), "both genomes")
    expect_error(getOrthoXML(), "valid HOG ID")
    expect_error(annotateSequence(), "sequence to query")
    expect_error(annotateSequence(c("A", "B")), "single sequence")
    expect_error(getHOG("HOG:1", level = 1), "taxonomic level")
    expect_error(getTaxonomy(newick = NA), "TRUE or FALSE")
    expect_length(calls$urls, 0L)
})
