#' Retrieve an OMA Group from the OMA Browser
#'
#' This function obtains an OMA Group from the OMA Browser database. An OMA
#' Group is defined to be a clique of proteins that are all orthologous to
#' each other, i.e. they are all related through speciation events only. An
#' OMA Group can thus by definition not contain any inparalogs. It is a very
#' stringent orthology grouping approach. OMA Groups are mostly useful to
#' infer phylogenetic species tree where they can be used as marker genes.
#'
#' Retrieving an OMA Group can be done using a group nr as id, its
#' fingerprint (a 7mer AA sequence which is unique to proteins in that
#' group), a member protein id or any sequence pattern that is unique to the
#' group.
#'
#' @param id An identifier for the group. See above for possible types of
#'     IDs.
#' @param attribute an extra attribute to be returned ('close_groups', the
#'     sorted list of closely related groups).
#' @return an \code{\link{omadb_obj}} with the group's attributes (the
#'     members are a data.frame), or a data.frame with the closely related
#'     groups.
#' @export
#' @examples
#' grp <- getOMAGroup(id = "YEAST58")
#' grp$fingerprint
#' head(grp$members[, c("omaid", "canonicalid")])
#' getOMAGroup(id = "ACCQRNQ")$group_nr
getOMAGroup <- function(id, attribute = NULL) {
    if (missing(id)) {
        stop("You must provide a valid group id.")
    }
    .checkScalarId(id)
    if (!is.null(attribute)) {
        attribute <- match.arg(attribute, "close_groups")
    }
    .toR(.omaGet(.apiUrl("group", id, attribute)))
}
