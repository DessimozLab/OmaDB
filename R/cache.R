## Two-level cache for API responses, similar to the one of the pyomadb
## package: an in-memory cache that is always active, and an optional
## persistent on-disk cache. As the content of the OMA Browser changes only
## with a new release, the persistent cache is cleared automatically when
## the OMA release or the API version reported by the server change.

.memCache <- function() {
    if (is.null(.omadb$memCache)) {
        size <- getOption("OmaDB.memCacheSize", 256 * 1024^2)
        .omadb$memCache <- cachem::cache_mem(max_size = size)
    }
    .omadb$memCache
}

#' Configure the caching of OMA Browser API responses
#'
#' All responses from the OMA Browser API are cached in memory for the
#' duration of the R session, so repeated queries (including accessing
#' attributes of objects that are loaded lazily) do not hit the server
#' again. Optionally, responses can also be stored in a persistent cache on
#' disk, which is reused across R sessions.
#'
#' The persistent cache is stored by default in the user cache directory
#' given by \code{tools::R_user_dir("OmaDB", which = "cache")}. As the data
#' in the OMA Browser only changes with a new release, the persistent cache
#' is cleared automatically if the OMA release or the API version of the
#' server differ from those of the cached data.
#'
#' The persistent cache can also be enabled by setting the option
#' \code{OmaDB.persistentCache} to \code{TRUE} before the first request,
#' e.g. in the \code{.Rprofile} file. The option \code{OmaDB.cacheDir} can
#' then be used to set its directory.
#'
#' @param persistent logical, whether responses should also be cached on
#'     disk.
#' @param dir directory of the persistent cache.
#' @param max_size maximum size of the persistent cache in bytes.
#' @return \code{setOmaCache} returns the directory of the persistent cache
#'     (or \code{NULL} if it is disabled) invisibly. \code{clearOmaCache}
#'     returns \code{NULL} invisibly.
#' @export
#' @examples
#' ## use a temporary directory for the persistent cache
#' setOmaCache(persistent = TRUE, dir = file.path(tempdir(), "omadb"))
#' clearOmaCache()
#' setOmaCache(persistent = FALSE)
setOmaCache <- function(persistent = TRUE,
        dir = tools::R_user_dir("OmaDB", which = "cache"),
        max_size = 1024^3) {
    stopifnot(is.logical(persistent), length(persistent) == 1L,
        !is.na(persistent))
    .omadb$diskChecked <- FALSE
    if (!persistent) {
        .omadb$diskCache <- NULL
        return(invisible(NULL))
    }
    stopifnot(is.character(dir), length(dir) == 1L)
    .omadb$diskCache <- cachem::cache_disk(dir = dir, max_size = max_size,
        evict = "lru")
    invisible(dir)
}

#' @rdname setOmaCache
#' @export
clearOmaCache <- function() {
    .memCache()$reset()
    disk <- .diskCache()
    if (!is.null(disk)) {
        disk$reset()
        .omadb$diskChecked <- FALSE
    }
    invisible(NULL)
}

.diskCache <- function() {
    if (is.null(.omadb$diskCache) &&
            isTRUE(getOption("OmaDB.persistentCache", FALSE))) {
        dir <- getOption("OmaDB.cacheDir",
            tools::R_user_dir("OmaDB", which = "cache"))
        setOmaCache(persistent = TRUE, dir = dir)
    }
    .omadb$diskCache
}

## Make sure the persistent cache belongs to the current OMA release.
##
## The release information is fetched uncached once per session (and after
## changing the cache settings) and compared with the release stored next to
## the cached data.
.checkDiskCache <- function(disk) {
    if (.omadb$diskChecked) {
        return(invisible(TRUE))
    }
    .omadb$diskChecked <- TRUE
    url <- .apiUrl("version")
    resp <- .perform(.buildRequest(url))
    current <- httr2::resp_body_json(resp)
    current <- paste(.omadb$apiUrl, current$oma_version,
        current$api_version)
    stored <- disk$get("omareleaseinfo", missing = NULL)
    if (!identical(stored, current)) {
        if (!is.null(stored)) {
            message("OMA release or API version changed; ",
                "clearing the persistent cache.")
        }
        disk$reset()
        disk$set("omareleaseinfo", current)
    }
    invisible(TRUE)
}

## Cache keys must consist of lowercase letters and digits only.
.cacheKey <- function(url, body = NULL) {
    rlang::hash(list(url, body, .omadb$apiVersion))
}

.cacheGet <- function(key) {
    mem <- .memCache()
    val <- mem$get(key, missing = NULL)
    if (!is.null(val)) {
        return(val)
    }
    disk <- .diskCache()
    if (!is.null(disk)) {
        .checkDiskCache(disk)
        val <- disk$get(key, missing = NULL)
        if (!is.null(val)) {
            mem$set(key, val)
        }
    }
    val
}

.cacheSet <- function(key, value) {
    .memCache()$set(key, value)
    disk <- .diskCache()
    if (!is.null(disk)) {
        .checkDiskCache(disk)
        disk$set(key, value)
    }
    invisible(value)
}
