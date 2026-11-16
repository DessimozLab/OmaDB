## Helpers to run tests against a fake OMA Browser API (no network access).

## Build an httr2 response for a fake API route.
fakeResponse <- function(url, body, status = 200L,
        type = "application/json", headers = list()) {
    httr2::response(status_code = status, url = url,
        headers = c(list(`Content-Type` = type), headers),
        body = charToRaw(body))
}

## Mock the OMA Browser API for the calling test.
##
## `routes` is a named list mapping urls to either a JSON string or a list
## with the elements body, status, type and headers. Returns an
## environment whose element `urls` records the requested urls.
local_fake_api <- function(routes, env = parent.frame()) {
    calls <- new.env()
    calls$urls <- character()
    mock <- function(req) {
        url <- req$url
        calls$urls <- c(calls$urls, url)
        route <- routes[[url]]
        if (is.null(route)) {
            return(fakeResponse(url, '{"detail": "not mocked"}',
                status = 404L))
        }
        if (is.character(route)) {
            route <- list(body = route)
        }
        fakeResponse(url, route$body,
            status = if (is.null(route$status)) 200L else route$status,
            type = if (is.null(route$type)) "application/json" else
                route$type,
            headers = if (is.null(route$headers)) list() else
                route$headers)
    }
    httr2::local_mocked_responses(mock, env = env)
    local_clean_state(env)
    calls
}

## Reset the API url and caches, and restore them at the end of the test.
local_clean_state <- function(env = parent.frame()) {
    oldUrl <- setAPI()
    setOmaCache(persistent = FALSE)
    clearOmaCache()
    withr::local_options(OmaDB.persistentCache = NULL, OmaDB.cacheDir = NULL,
        .local_envir = env)
    withr::defer({
        setAPI(oldUrl)
        setOmaCache(persistent = FALSE)
        clearOmaCache()
    }, envir = env)
    invisible(NULL)
}

api <- function(path) {
    paste0("https://omabrowser.org/api/", path)
}
