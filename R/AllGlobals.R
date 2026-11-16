## Package-wide state: the API configuration (see client.R) and the
## response caches (see cache.R).
.omadb <- new.env(parent = emptyenv())
.omadb$defaultUrl <- "https://omabrowser.org/api/"
.omadb$apiUrl <- .omadb$defaultUrl
.omadb$apiVersion <- "1.11"
.omadb$perPage <- 5000L
.omadb$memCache <- NULL
.omadb$diskCache <- NULL
.omadb$diskChecked <- FALSE
