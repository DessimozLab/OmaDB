## Integration tests against recorded responses of the OMA Browser API,
## see helper-recorded-api.R.

test_that("getVersion returns the release information", {
    with_recorded_api({
        v <- getVersion()
        expect_s3_class(v, "omadb_obj")
        expect_match(v$oma_version, "^All\\.")
        expect_type(v$api_version, "character")
    })
})

test_that("getProtein returns a protein object", {
    with_recorded_api({
        prot <- getProtein("YEAST00001")
        expect_s3_class(prot, "omadb_obj")
        expect_identical(prot$omaid, "YEAST00001")
        expect_identical(prot$species$code, "YEAST")
        expect_s4_class(prot$sequence, "AAString")
        expect_s4_class(prot$cdna, "DNAString")
        expect_s4_class(prot$locus, "GRanges")
        expect_identical(length(prot$sequence), prot$sequence_length)
        orth <- prot$orthologs
        expect_s3_class(orth, "data.frame")
        expect_true(all(c("omaid", "rel_type", "species.code") %in%
            names(orth)))
    })
})

test_that("getProtein retrieves attributes of proteins", {
    with_recorded_api({
        orth <- getProtein("YEAST00001", attribute = "orthologs",
            rel_type = "1:1")
        expect_s3_class(orth, "data.frame")
        expect_true(all(orth$rel_type == "1:1"))

        xref <- getProtein("YEAST00001", attribute = "xref")
        expect_true("P39728" %in% xref$xref)

        dom <- getProtein("YEAST00001", attribute = "domains")
        expect_s3_class(dom, "omadb_obj")
        expect_s3_class(dom$regions, "data.frame")

        iso <- getProtein("YEAST00001", attribute = "isoforms")
        expect_s3_class(iso, "data.frame")

        locus <- getProtein("YEAST00001", attribute = "locus")
        expect_s4_class(locus, "GRanges")
    })
})

test_that("getProtein retrieves several proteins at once", {
    with_recorded_api({
        ids <- c("YEAST00001", "YEAST00002")
        prots <- getProtein(ids)
        expect_type(prots, "list")
        expect_named(prots, ids)
        expect_identical(vapply(prots, function(p) p$omaid, ""),
            stats::setNames(ids, ids))

        go <- getProtein(ids, attribute = "gene_ontology")
        expect_s3_class(go, "data.frame")
        expect_setequal(unique(go$query_id), ids)
        expect_true(all(grepl("^GO:", go$GO_term)))

        loci <- getLocus(ids)
        expect_s4_class(loci, "GRanges")
        expect_identical(names(loci), ids)
    })
})

test_that("getHOG returns HOGs, their members and annotations", {
    with_recorded_api({
        hog <- getHOG("YEAST00001")
        expect_s3_class(hog, "data.frame")
        expect_identical(nrow(hog), 1L)
        expect_match(hog$hog_id, "^HOG:")

        root <- getHOG("YEAST00001", level = "root")
        expect_match(root$hog_id, "^HOG:F[0-9]+$")

        members <- getHOG("YEAST00001", level = "Saccharomycetaceae",
            members = TRUE)
        expect_s3_class(members, "data.frame")
        expect_true("YEAST00001" %in% members$omaid)
        expect_identical(getHOG("YEAST00001", level = "Saccharomycetaceae",
            attribute = "members"), members)

        go <- getHOG("YEAST00001", attribute = "gene_ontology")
        expect_s3_class(go, "data.frame")
        expect_true(all(c("GO_term", "score") %in% names(go)))

        sim <- getHOG("YEAST00001", attribute = "similar_profile_hogs",
            max_results = 3)
        expect_s3_class(sim, "data.frame")
        expect_lte(nrow(sim), 3L)
    })
})

test_that("getOMAGroup returns a group with its members", {
    with_recorded_api({
        grp <- getOMAGroup("YEAST00058")
        expect_s3_class(grp, "omadb_obj")
        expect_s3_class(grp$members, "data.frame")
        expect_true("YEAST00058" %in% grp$members$omaid)
        loci <- getLocus(grp$members)
        expect_s4_class(loci, "GRanges")
        expect_identical(length(loci), nrow(grp$members))
    })
})

test_that("getGenome returns genomes and their proteins", {
    with_recorded_api({
        genome <- getGenome("ASHGO")
        expect_s3_class(genome, "omadb_obj")
        expect_identical(genome$code, "ASHGO")
        expect_s3_class(genome$chromosomes, "data.frame")
    })
})

test_that("getGenomePairs and getSharedAncestry compare genomes", {
    with_recorded_api({
        pairs <- getGenomePairs("YEAST", "ASHGO", chr1 = "I",
            rel_type = "1:1")
        expect_s3_class(pairs, "data.frame")
        expect_gt(nrow(pairs), 0L)
        expect_true(all(pairs$entry_1.chromosome == "I"))
        expect_true(all(pairs$rel_type == "1:1"))
        expect_true(all(startsWith(pairs$entry_1.omaid, "YEAST")))
        expect_true(all(startsWith(pairs$entry_2.omaid, "ASHGO")))

        shared <- getSharedAncestry("YEAST", "ASHGO")
        expect_gt(shared$fraction, 0)
        expect_lte(shared$fraction, 1)
        expect_s3_class(shared$details, "data.frame")
    })
})

test_that("getTaxonomy and getTree return the taxonomy", {
    with_recorded_api({
        tax <- getTaxonomy(members = c("YEAST", "ASHGO"))
        expect_type(tax$newick, "character")
        tree <- getTree(tax)
        expect_s3_class(tree, "phylo")
        expect_length(tree$tip.label, 2L)
        expect_true(any(grepl("(strain", tree$tip.label, fixed = TRUE)))

        dict <- getTaxonomy(members = c("YEAST", "ASHGO"), newick = FALSE)
        expect_s3_class(dict, "omadb_obj")
        expect_s3_class(dict$children, "data.frame")
    })
})

test_that("sequences can be mapped, annotated and searched", {
    seq <- paste0("MNDPSLLGYPNVGPQQQQQQQQQQHAGLLGKGTPNALQQQLHMNQLTGIPPPGL",
        "MNNSDVHTSSNNNSRQLLDQLANGNANMLNMNMDNNNNNNNNNNNNNNNGGGSGVMMNASTA")
    with_recorded_api({
        res <- mapSequence(query = Biostrings::AAString(seq))
        expect_s3_class(res, "omadb_obj")
        expect_s3_class(res$targets, "data.frame")
        expect_gte(nrow(res$targets), 1L)

        go <- annotateSequence(seq)
        expect_s3_class(go, "data.frame")
        expect_true("GO_ID" %in% names(go))

        hits <- searchProtein("P53_RAT")
        expect_s3_class(hits, "data.frame")
        expect_true("P53_RAT" %in% hits$xref)
    })
})

test_that("getSynteny returns a graph as nodes and links", {
    with_recorded_api({
        syn <- getSynteny(level = "YEAST", id = "YEAST00001", context = 1)
        expect_named(syn, c("nodes", "links"))
        expect_true("YEAST00001" %in% syn$nodes$id)
        expect_identical(names(syn$links)[1:2], c("source", "target"))
        expect_identical(names(syn$nodes)[1], "id")
    })
})

test_that("getOrthoXML returns an OrthoXML document", {
    with_recorded_api({
        xml <- getOrthoXML("YEAST00001")
        expect_type(xml, "character")
        expect_match(xml, "<orthoXML")
    })
})

test_that("formatTopGO converts GO annotations", {
    with_recorded_api({
        prots <- getProtein(c("YEAST00001", "YEAST00002"))
        g2go <- formatTopGO(prots, format = "geneID2GO")
        expect_named(g2go, c("YEAST00001", "YEAST00002"))
        expect_true(all(grepl("^GO:", unlist(g2go))))

        go <- getProtein(c("YEAST00001", "YEAST00002"),
            attribute = "gene_ontology")
        expect_identical(formatTopGO(go, format = "geneID2GO"), g2go)
        go2g <- formatTopGO(go, format = "GO2geneID")
        expect_true(all(grepl("^GO:", names(go2g))))
    })
})
