## Low-level access to the OMA Browser REST API.
##
## All exported functions build an API path and call .omaGet() / .omaPost(),
## which take care of caching (see cache.R), retries, pagination and error
## reporting. Responses are parsed with jsonlite into plain lists and then
## converted into R objects by the helpers in convert.R.

.userAgent <- function() {
    paste0("r-omadb/", utils::packageVersion("OmaDB"))
}

#' Set the url to the OMA Browser API
#'
#' Sets the base url of the OMA Browser REST API used by all functions of
#' the package. If no url is given, the default OMA Browser API
#' (\url{https://omabrowser.org/api/}) is used. This can be used to access
#' other instances of the OMA Browser, e.g. the Corona OMA Browser at
#' \url{https://corona.omabrowser.org/api/}.
#'
#' Changing the API url does not invalidate cached responses, as the url is
#' part of the cache key.
#'
#' @param url Base url to the API
#' @return The previously used API url, invisibly.
#' @export
#' @examples
#' old <- setAPI("https://corona.omabrowser.org/api/")
#' setAPI(old)
setAPI <- function(url) {
    if (missing(url) || is.null(url)) {
        url <- .omadb$defaultUrl
    }
    stopifnot(is.character(url), length(url) == 1L, !is.na(url))
    if (!endsWith(url, "/")) {
        url <- paste0(url, "/")
    }
    old <- .omadb$apiUrl
    .omadb$apiUrl <- url
    invisible(old)
}

## Build the url of an API endpoint.
##
## `...` are path components (NULL components are dropped, all others are
## url-encoded), `query` is a named list of query parameters (NULL values
## are dropped).
.apiUrl <- function(..., query = list()) {
    parts <- Filter(Negate(is.null), list(...))
    parts <- vapply(parts, function(p) {
        utils::URLencode(as.character(p), reserved = TRUE)
    }, character(1))
    url <- paste0(.omadb$apiUrl, paste0(parts, "/", collapse = ""))
    .addQuery(url, query)
}

.addQuery <- function(url, query) {
    query <- Filter(Negate(is.null), query)
    if (length(query) == 0L) {
        return(url)
    }
    vals <- vapply(query, function(v) {
        if (is.logical(v)) {
            v <- tolower(as.character(v))
        }
        utils::URLencode(paste(as.character(v), collapse = ","),
            reserved = TRUE)
    }, character(1))
    sep <- if (grepl("?", url, fixed = TRUE)) "&" else "?"
    paste0(url, sep, paste(names(query), vals, sep = "=", collapse = "&"))
}

## Perform a single HTTP request and return the response body as text.
##
## Successful responses are cached. The return value is a list with the
## body text and the url of the next page (or NULL) for paginated
## endpoints.
.fetch <- function(url, body = NULL, accept = NULL) {
    key <- .cacheKey(url, body)
    hit <- .cacheGet(key)
    if (!is.null(hit)) {
        return(hit)
    }
    resp <- .perform(.buildRequest(url, body, accept))
    res <- list(text = httr2::resp_body_string(resp),
        nextUrl = .nextLink(resp))
    .cacheSet(key, res)
    res
}

## All requests go through this function (allows mocking in unit tests).
.perform <- function(req) {
    httr2::req_perform(req)
}

.buildRequest <- function(url, body = NULL, accept = NULL) {
    if (is.null(accept)) {
        accept <- paste0("application/json; version=", .omadb$apiVersion)
    }
    req <- httr2::request(url)
    req <- httr2::req_user_agent(req, .userAgent())
    req <- httr2::req_headers(req, Accept = accept)
    req <- httr2::req_timeout(req, getOption("timeout", 60))
    req <- httr2::req_retry(req, max_tries = 3L,
        is_transient = .isTransient, backoff = function(i) 2^i)
    req <- httr2::req_error(req, body = .errorBody)
    if (!is.null(body)) {
        req <- httr2::req_body_raw(req, body, type = "application/json")
    }
    req
}

.isTransient <- function(resp) {
    httr2::resp_status(resp) %in% c(429L, 502L, 503L, 504L)
}

## Additional lines of the error message for failed requests.
.errorBody <- function(resp) {
    url <- httr2::resp_url(resp)
    type <- httr2::resp_content_type(resp)
    detail <- NULL
    if (identical(type, "application/json")) {
        content <- tryCatch(httr2::resp_body_json(resp),
            error = function(e) NULL)
        if (.isRecord(content) && !is.null(content$detail)) {
            detail <- as.character(content$detail)
        } else if (.isScalarList(content)) {
            detail <- as.character(unlist(content))
        }
    } else if (identical(httr2::resp_status(resp), 403L) &&
            grepl("text/html", type, fixed = TRUE)) {
        detail <- paste("The request was blocked by the web application",
            "firewall of the OMA Browser.")
    }
    c(detail, paste("Request url:", url))
}

.nextLink <- function(resp) {
    if (!httr2::resp_header_exists(resp, "link")) {
        return(NULL)
    }
    httr2::resp_link_url(resp, "next")
}

.parseJSON <- function(text) {
    jsonlite::parse_json(text, simplifyVector = FALSE)
}

## GET an endpoint and return the parsed JSON.
##
## If the endpoint is paginated, all pages are retrieved and concatenated.
.omaGet <- function(url, paginate = TRUE) {
    res <- .fetch(url)
    data <- .parseJSON(res$text)
    while (paginate && !is.null(res$nextUrl)) {
        res <- .fetch(res$nextUrl)
        data <- c(data, .parseJSON(res$text))
    }
    data
}

## GET a paginated list endpoint using the package's page size.
.omaGetList <- function(url) {
    if (!grepl("[?&]per_page=", url)) {
        perPage <- getOption("OmaDB.perPage", .omadb$perPage)
        url <- .addQuery(url, list(per_page = perPage))
    }
    .omaGet(url)
}

## POST a json body to an endpoint and return the parsed JSON.
.omaPost <- function(url, body) {
    .parseJSON(.fetch(url, body = body)$text)
}

## Is `x` a url that points to the configured OMA Browser API?
.isApiUrl <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) &&
        grepl("^https?://", x) && grepl("/api/", x, fixed = TRUE)
}
