## template-depr_fun.r
#' @name <%= fun %>-deprecated
#' @usage <%= gsub("\n", "\n#' ", gsub("\u{A0}", " ", roxygen2:::wrapString(roxygen2:::function_usage(fun, formals(fun))))) %>
#' @seealso \code{\link{OmaDB-deprecated}}
#' @keywords internal
