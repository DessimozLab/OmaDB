
#' @importFrom utils URLencode packageVersion
#' @import httr
#' @import plyr
#' @import ape
#' @importFrom  Biostrings AAString
#' @importFrom Biostrings DNAString
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom IRanges IRanges


pkg.env <- new.env()
pkg.env$API_default_url <- "https://omabrowser.org/api/"
pkg.env$API_url <- pkg.env$API_default_url
pkg.env$API_version <- "1.6"
pkg.env$user_agent <- paste0("r-omadb/", packageVersion("OmaDB"))


urlGenerator <- function(endpoint, id = NULL, detail = NULL, ...) {
    url_prefix <- paste0(pkg.env$API_url, tolower(endpoint), "/")
    if (!is.null(id)) {
        id <- paste0(utils::URLencode(as.character(id)), "/")
    }
    if (!is.null(detail)) {
        detail <- paste0(as.character(detail), "/")
    }
    base_url <- paste0(url_prefix, id, detail)

    qargs <- unlist(list(...))
    if (length(qargs) == 0) {
        return(base_url)
    }

    # remove all qargs that are NULL
    qargs[sapply(qargs, is.null)] <- NULL
    if (length(qargs) == 0) {
        return(base_url)
    }

    qarg_vals <- lapply(lapply(qargs, as.character), URLencode)
    qargs <- mapply(paste, names(qargs), qarg_vals, sep = "=", USE.NAMES = F)
    query_param <- paste(qargs, collapse = "&")
    final_url <- paste(base_url, query_param, sep = "?")
    return(final_url)
}


simpleRequest <- function(url, body = NULL) {
    response <- load_url(url, body = body)
    if (!is.null(response)) {
        content_list <- httr::content(response, as = "parsed")
        column_names <- names(content_list)
        return(objectFactory(column_names, content_list))
    }
    return(NULL)
}

load_url <- function(url, body = NULL) {
    count <- 0
    accept_header = paste0("application/json; version=", pkg.env$API_version)
    while (count < 3) {
        response <- tryCatch(if (!is.null(body)) {
            httr::POST(url, body = body, encode = "raw", httr::accept(accept_header),
                httr::content_type("application/json"), httr::user_agent(pkg.env$user_agent))
        } else {
            httr::GET(url, httr::accept(accept_header), httr::content_type("application/json"),
                httr::user_agent(pkg.env$user_agent))
        }, error = function(cond) {
            Sys.sleep(0.5)
            return(NULL)
        }, warning = function(cond) {
            message(cond)
            return(NULL)
        })

        if (is.null(response)) {
            count <- count + 1
            next
        }

        if (response$status_code < 500) {
            break
        }
        message(sprintf("Request failed with server error %d. retry in %d secs",
            response$status_code, 2^count))
        Sys.sleep(2^count)
        count <- count + 1
    }

    if (is.null(response)) {
        message(sprintf("Cannot resolve '%s'. Please connect to the internet", url))
        return(NULL)
    }

    tryCatch({
        httr::stop_for_status(response)
        # request worked out
        return(response)
    }, http_400 = function(cond) {
        reason <- httr::content(response, as = "parsed")
        if (!is.null(reason)) {
            cond$message <- paste0(cond$message, " Reason: ", reason$detail)
        }
        message(cond)
        return(NULL)
    }, error = function(cond) {
        cond$message <- paste0("Request to ", url, " failed:\n", cond$message, "\n")
        message(cond)
        return(NULL)
    }, warning = function(cond) {
        warning(cond)
        return(NULL)
    })
}


objectFactory <- function(column_names, content_list) {
    list_of_variables <- lapply(column_names, FUN = function(name) {
        content <- content_list[[name]]
        if (is.null(content)) {
            content == ""
        }
        if (class(content) == "list" && name %in% c("children_hogs", "alternative_levels",
            "parent_hogs")) {
            # delete these fields
            return(NULL)
        }

        if (class(content) == "list" && length(content) != 0 && name != "locus" &&
            name != "chromosomes") {
            if (is.null(names(content))) {
                formatData(content)  # if list, make into DF.
            }
        } else if (name == "chromosomes") {
            formatData(content)
        } else if (name == "locus") {
            GenomicRanges::GRanges(seqnames = content_list[["omaid"]], ranges = IRanges::IRanges(content$start,
                content$end), strand = content$strand)
        } else if (name == "sequence") {
            Biostrings::AAString(content)
        } else if (name == "cdna") {
            content <- gsub("X", "N", content)
            Biostrings::DNAString(content)
        } else {
            content
        }
    })

    names(list_of_variables) <- column_names
    list_of_variables[sapply(list_of_variables, is.null)] <- NULL
    value <- list_of_variables
    class(value) <- "omadb_obj"
    return(value)
}

largeRequestFactory <- function(url, n, per_page, body = NULL) {
    if (!is.null(body)) {
        response_list <- lapply(body, FUN = function(x) {
            load_url(url, body = x)
        })
    } else {
        n_requests <- ceiling(as.numeric(n) / per_page)
        prefix <- if (substr(url, nchar(url), nchar(url)) == "/")
            "?per_page=" else "&per_page="
        url_list <- list()
        for (i in seq_along(1:n_requests)) {
            url_list[[i]] <- paste0(url, prefix, per_page, "&page=", i)
        }
        response_list <- lapply(url_list, FUN = function(x) {
            load_url(x, body)
        })
    }

    json_list <- lapply(response_list, FUN = function(x) {
        httr::content(x, as = "parsed")
    })
    content_list <- do.call("c", json_list)
    return(content_list)
}

extractdata <- function(content_list) {
    column_names <- names(content_list)
    if (is.null(column_names)) {
        # this is due to /hog/ formatting in the API
        if (length(content_list) == 0) {
            return(content_list)
        }

        if (length(content_list) == 1) {
            column_names <- names(content_list[[1]])
            content_list <- content_list[[1]]
            return(objectFactory(column_names, content_list))

        } else if (length(content_list) > 1) {
            return(formatData(content_list))
        }
    } else {
        return(objectFactory(column_names, content_list))
    }
}



requestFactory <- function(url, body = NULL, per_page = 5000, page = NULL) {
    if (!is.null(body)) {
        return(largeRequestFactory(url, body = body))
    }
    
    # sep for per_page query params is either ? or & depending if 
    # no query param so far or not
    sep <- if (substr(url, nchar(url), nchar(url)) == "/")
        "?" else "&"
    qq <- paste0("per_page=", per_page, "&page=", if (is.null(page))
        1 else page)
    first_url <- paste(url, qq, sep = sep)

    response <- load_url(first_url, body = body)
    if (is.null(response)) {
        return(NULL)
    }

    content_list <- httr::content(response, as = "parsed")
    n_items <- httr::headers(response)[["x-total-count"]]
    if (is.null(page) && !is.null(n_items)) {
        # if we need all data and data is pageinated, check if we already have all data.
        # if not, retrieve in parallel
        nr_elem <- length(content_list)
        n_items <- as.numeric(n_items)
        if (nr_elem < n_items) {
            content_list <- largeRequestFactory(url, n_items, per_page)
        }
    }
    return(extractdata(content_list))
}



formatData <- function(data) {
    # flatten chromosomes in genome so that GRanges object can form
    if ("entry_ranges" %in% names(data[[1]])) {
        for (i in seq_along(data)) {
            data[[i]][["entry_ranges.1"]] <- data[[i]]["entry_ranges"][[1]][[1]][[1]]
            data[[i]][["entry_ranges.2"]] <- data[[i]]["entry_ranges"][[1]][[1]][[2]]
            data[[i]][["entry_ranges"]] <- NULL
        }
    } else if ("sequence" %in% names(data[[1]]) || "alternative_levels" %in% names(data[[1]])) {
        # to ensure bulk retrieve nor hog-info does not form a df
        for (i in seq_along(data)) {
            data[[i]] <- objectFactory(names(data[[i]]), data[[i]])
        }
        return(data)
    }

    dfs <- lapply(data, data.frame, stringsAsFactors = FALSE)
    data <- plyr::rbind.fill(dfs)
    return(data)
}

#' Resolve URLs automatically when accessed
#'
#' The function to obtain further information from a given url.
#'
#' @param x object
#' @param name attribute
#' @return API response behind the URL
#' @export

"$.omadb_obj" <- function(x, name) {
    if (is.character(x[[name]]) && grepl("https://", x[[name]])) {
        value <- requestFactory(x[[name]])
        obj_name <- deparse(substitute(x))
        set_new_val(x, name) <- value
        assign(obj_name, x, envir = .GlobalEnv)
    }
    return(x[[name]])
}

"set_new_val<-" <- function(x, name, value) {
    x[[name]] <- value
    x
}


#' Set the url to the OMA Browser API
#'
#' Function to set the base url to the OMA Browser API. If no url is
#' specified, the default OMA Browser API url is used.
#'
#' @param url Base url to the API
#' @export

setAPI <- function(url) {
    if (missing(url)) {
        url <- pkg.env$API_default_url
    }
    if (substr(url, nchar(url), nchar(url)) != "/") {
        url <- paste0(url, "/")
    }
    pkg.env$API_url <- url
}


#' Load data for a given url from the OMA Browser API.
#'
#' This function is usualy not needed by users. In most circumstances
#' an attribute containing a URL is automatically loaded when accessed.
#' However, in case the data is transformed into a dataframe, this will
#' no longer be true, in which case one can access the data behind this
#' attribute using this function.
#'
#' @param url The url of interest
#' @return a data.frame containing the information behind an URL
#' @export
#' @examples
#' resolveURL('http://omabrowser.org/api/protein/YEAST58/gene_ontology/')
resolveURL <- function(url) {
    return(requestFactory(url))
}
