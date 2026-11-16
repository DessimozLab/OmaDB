test_that("urls are built and encoded correctly", {
    local_clean_state()
    expect_identical(.apiUrl("protein", "P53_RAT"), api("protein/P53_RAT/"))
    expect_identical(.apiUrl("hog", "HOG:F0411479", "members"),
        api("hog/HOG%3AF0411479/members/"))
    expect_identical(.apiUrl("hog", NULL, query = list(level = "Homo sapiens",
        compare_with = NULL)), api("hog/?level=Homo%20sapiens"))
    expect_identical(.apiUrl("sequence", query = list(query = "MKV",
        full_length = FALSE)), api("sequence/?query=MKV&full_length=false"))
    expect_identical(.addQuery(api("x/?a=1"), list(b = c("A", "B"))),
        api("x/?a=1&b=A%2CB"))
})

test_that("setAPI changes the base url", {
    local_clean_state()
    old <- setAPI("https://corona.omabrowser.org/api")
    expect_identical(old, "https://omabrowser.org/api/")
    expect_identical(.apiUrl("version"),
        "https://corona.omabrowser.org/api/version/")
    setAPI()
    expect_identical(.apiUrl("version"), api("version/"))
})

test_that("paginated endpoints are loaded completely", {
    withr::local_options(OmaDB.perPage = 2L)
    calls <- local_fake_api(list(
        "https://omabrowser.org/api/genome/?per_page=2" = list(
            body = '[{"code": "A"}, {"code": "B"}]',
            headers = list(Link = paste0(
                '<https://omabrowser.org/api/genome/?page=2&per_page=2>; ',
                'rel="next"'))),
        "https://omabrowser.org/api/genome/?page=2&per_page=2" =
            '[{"code": "C"}]'))
    genomes <- getGenome()
    expect_s3_class(genomes, "data.frame")
    expect_identical(genomes$code, c("A", "B", "C"))
    expect_length(calls$urls, 2L)
})

test_that("HTTP errors are reported with the server's explanation", {
    local_fake_api(list(
        "https://omabrowser.org/api/protein/NOTANID/" = list(status = 404L,
            body = '{"detail": "requested id \'NOTANID\' is unknown"}'),
        "https://omabrowser.org/api/xref/?search=MAL" = list(status = 400L,
            body = '["Query too unspecific. Matches >1000 elements"]'),
        "https://omabrowser.org/api/group/1/close_groups/" = list(
            status = 403L, type = "text/html",
            body = "<html>Just a moment...</html>")))
    expect_error(getProtein("NOTANID"), "is unknown",
        class = "httr2_http_404")
    expect_error(getProtein("NOTANID"), "protein/NOTANID")
    expect_error(searchProtein("MAL"), "too unspecific")
    expect_error(getOMAGroup(1, attribute = "close_groups"), "firewall",
        class = "httr2_http_403")
})

test_that("invalid arguments are rejected before querying the server", {
    calls <- local_fake_api(list())
    expect_error(getProtein(), "at least one protein id")
    expect_error(getProtein("P1", attribute = "foo"), "should be one of")
    expect_error(getProtein("P1", attribute = "orthologs", rel_type = "x"),
        "rel_type must be one of")
    expect_error(getHOG(), "valid HOG ID")
    expect_error(getHOG("HOG:1", members = "yes"), "TRUE or FALSE")
    expect_error(getGenome(attribute = "proteins"), "genome id")
    expect_error(getGenomePairs("YEAST"), "both genomes")
    expect_error(getTaxonomy(root = "A", members = "B"), "Only one of")
    expect_error(mapSequence(), "sequence to query")
    expect_error(mapSequence("MKV", search = "fuzzy"), "should be one of")
    expect_error(searchProtein("AB"), "at least 3 characters")
    expect_error(getSynteny(), "taxonomic level")
    expect_error(getHOGs(), "taxonomic level")
    expect_length(calls$urls, 0L)
})

test_that("defunct functions raise an error", {
    expect_error(getXref("MAL"), class = "defunctError")
    expect_error(getData("protein", "P1"), "getProtein")
    expect_error(getAnnotation("MKV"), class = "defunctError")
    expect_error(getGenomeAlignment("A", "B"), class = "defunctError")
})
