versionRoute <- function(release) {
    sprintf('{"oma_version": "%s", "api_version": "1.11"}', release)
}

test_that("responses are cached in memory", {
    calls <- local_fake_api(list(
        "https://omabrowser.org/api/protein/P1/" = '{"omaid": "P1"}'))
    expect_identical(getProtein("P1")$omaid, "P1")
    expect_identical(getProtein("P1")$omaid, "P1")
    expect_length(calls$urls, 1L)

    clearOmaCache()
    getProtein("P1")
    expect_length(calls$urls, 2L)
})

test_that("lazy attributes are loaded once and the object is unchanged", {
    calls <- local_fake_api(list(
        "https://omabrowser.org/api/protein/P1/" = paste0('{"omaid": "P1",',
            '"orthologs": "https://omabrowser.org/api/protein/P1/orthologs/"}'),
        "https://omabrowser.org/api/protein/P1/orthologs/" =
            '[{"omaid": "P2", "rel_type": "1:1"}]'))
    prot <- getProtein("P1")
    orth <- prot$orthologs
    expect_s3_class(orth, "data.frame")
    expect_identical(orth$omaid, "P2")
    expect_identical(getAttribute(prot, "orthologs"), orth)
    expect_true(.isApiUrl(prot[["orthologs"]]))
    expect_length(calls$urls, 2L)
})

test_that("the persistent cache survives clearing the memory cache", {
    dir <- withr::local_tempdir()
    calls <- local_fake_api(list(
        "https://omabrowser.org/api/version/" = versionRoute("All.Jan2026"),
        "https://omabrowser.org/api/protein/P1/" = '{"omaid": "P1"}'))
    setOmaCache(persistent = TRUE, dir = dir)
    getProtein("P1")
    expect_identical(calls$urls, c(api("version/"), api("protein/P1/")))

    ## a new session: empty memory cache, release is checked again
    .memCache()$reset()
    setOmaCache(persistent = TRUE, dir = dir)
    getProtein("P1")
    expect_identical(calls$urls[3], api("version/"))
    expect_length(calls$urls, 3L)
})

test_that("the persistent cache is cleared for a new OMA release", {
    dir <- withr::local_tempdir()
    release <- "All.Jan2026"
    calls <- new.env()
    calls$urls <- character()
    local_clean_state()
    httr2::local_mocked_responses(function(req) {
        calls$urls <- c(calls$urls, req$url)
        body <- if (endsWith(req$url, "version/")) versionRoute(release)
            else '{"omaid": "P1"}'
        fakeResponse(req$url, body)
    })
    setOmaCache(persistent = TRUE, dir = dir)
    getProtein("P1")
    expect_length(calls$urls, 2L)

    release <- "All.Jul2026"
    .memCache()$reset()
    setOmaCache(persistent = TRUE, dir = dir)
    expect_message(getProtein("P1"), "clearing the persistent cache")
    expect_identical(calls$urls[3:4], c(api("version/"), api("protein/P1/")))
})

test_that("the persistent cache can be enabled with an option", {
    dir <- withr::local_tempdir()
    local_fake_api(list(
        "https://omabrowser.org/api/version/" = versionRoute("All.Jan2026"),
        "https://omabrowser.org/api/protein/P1/" = '{"omaid": "P1"}'))
    withr::local_options(OmaDB.persistentCache = TRUE, OmaDB.cacheDir = dir)
    getProtein("P1")
    expect_false(is.null(.omadb$diskCache))
    expect_true(length(list.files(dir)) > 0L)
})

test_that("POST requests are cached by their body", {
    calls <- local_fake_api(list())
    httr2::local_mocked_responses(function(req) {
        calls$urls <- c(calls$urls, req$url)
        body <- req$body$data
        if (is.raw(body)) {
            body <- rawToChar(body)
        }
        ids <- jsonlite::fromJSON(body)$ids
        body <- jsonlite::toJSON(lapply(ids, function(id) {
            list(query_id = id, target = if (id != "NO") list(omaid = id))
        }), auto_unbox = TRUE, null = "null")
        fakeResponse(req$url, as.character(body))
    })
    res <- getProtein(c("P1", "NO", "P2"))
    expect_named(res, c("P1", "NO", "P2"))
    expect_null(res$NO)
    expect_identical(res$P2$omaid, "P2")
    getProtein(c("P1", "NO", "P2"))
    expect_length(calls$urls, 1L)
    getProtein(c("P1", "P2"))
    expect_length(calls$urls, 2L)
})
