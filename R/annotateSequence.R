#' Map GO annotation to a sequence that is not available in the OMA Browser
#'
#' This function obtain Gene Ontology annotation for a given sequence that
#' does not need to exist in the OMA Browser so far. The query sequence will
#' analysed and a fast homology detection approach based on kmers will be used
#' to detect the closest sequences in OMA. GO annotations for these top hits
#' will be used to annotated the query sequence.
#'
#' @param query the sequence to be annotated, it can be either a string or an
#'        AAString object from the Biostrings package
#' @return a data.frame containing the the GO annotation information of the
#'        most similar protein to the query sequence
#' @export
#' @examples
#' annotateSequence(query='MNDPSLLGYPNVGPQQQQQQQQQQHAGLLGKGTPNALQQQLHMNQLTGIPPPGLMNNSDVHTSSNNNSRQLLDQLANGNANMLNMNMDNNNNNNNNNNNNNNNGGGSGVMMNASTAAVNSIGMVPTVGTPVNINVNASNPLLHPHLDDPSLLNNPIWKLQLHLAAVSAQSLGQPNIYARQNAMKKYLATQQAQQAQQQAQQQAQQQVPGPFGPGPQAAPPALQPTDFQQSHIAEASKSLVDCTKQALMEMADTLTDSKTAKKQQPTGDSTPSGTATNSAVSTPLTPKIELFANGKDEANQALLQHKKLSQYSIDEDDDIENRMVMPKDSKYDDQLWHALDLSNLQIFNISANIFKYDFLTRLYLNGNSLTELPAEIKNLSNLRVLDLSHNRLTSLPAELGSCFQLKYFYFFDNMVTTLPWEFGNLCNLQFLGVEGNPLEKQFLKILTEKSVTGLIFYLRDNRPEIPLPHERRFIEINTDGEPQREYDSLQQSTEHLATDLAKRTFTVLSYNTLCQHYATPKMYRYTPSWALSWDYRRNKLKEQILSYDSDLLCLQEVESKTFEEYWVPLLDKHGYTGIFHAKARAKTMHSKDSKKVDGCCIFFKRDQFKLITKDAMDFSGAWMKHKKFQRTEDYLNRAMNKDNVALFLKLQHIPSGDTIWAVTTHLHWDPKFNDVKTFQVGVLLDHLETLLKEETSHNFRQDIKKFPVLICGDFNSYINSAVYELINTGRVQIHQEGNGRDFGYMSEKNFSHNLALKSSYNCIGELPFTNFTPSFTDVIDYIWFSTHALRVRGLLGEVDPEYVSKFIGFPNDKFPSDHIPLLARFEFMKTNTGSKKV')

annotateSequence <- function(query) {
    if (missing(query)) {
        stop("You must provide a sequence to query.")
    }
    if (class(query) == "AAString") {
        query <- as.character(query)
    }

    url <- urlGenerator(endpoint = "function", query = query)
    return(requestFactory(url))
}


