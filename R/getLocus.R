#' Get loci for a given list of proteins
#'
#' Function to obtain loci in genomic range format for a given list of proteins
#'
#' @param proteins the dataframe or a list of dataframes containing the protein data of interest. this can either be the members df or a list of protein ids.
#' @return genomic range object from the GenomicRanges package in Bioconductor
#' @export
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @examples
#' loci = getLocus(proteins = getOMAGroup('YEAST58')['members'])

getLocus <- function(proteins) {

    # if a list of dataframes is given.
    if (class(proteins) == "list") {
        if (class(proteins[[1]]) == "data.frame") {

            df <- plyr::rbind.fill(proteins)

        } else {

            df <- getProtein(id = proteins, attribute = "locus")

        }

    }

    df$locus.strand[df$locus.strand == "1"] <- "+"
    df$locus.strand[df$locus.strand == "-1"] <- "-"

    grange <- GenomicRanges::makeGRangesFromDataFrame(df, start.field = "start",
        end.field = "end", seqnames.field = "omaid")

    return(grange)

}

