## Integration tests replay responses of the OMA Browser API that were
## recorded with httptest2 in tests/testthat/omadb/. To re-record them
## (e.g. after a new OMA release), run the tests with network access and
## the environment variable OMADB_RECORD_API set to "true".

with_recorded_api <- function(expr) {
    local_clean_state()
    oldPaths <- httptest2::.mockPaths()
    options(httptest2.mock.paths = testthat::test_path("omadb"))
    on.exit(options(httptest2.mock.paths = oldPaths), add = TRUE)
    if (identical(Sys.getenv("OMADB_RECORD_API"), "true")) {
        httptest2::capture_requests(expr, simplify = FALSE)
    } else {
        httptest2::with_mock_api(expr)
    }
}
