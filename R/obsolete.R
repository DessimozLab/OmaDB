getXref <- function(pattern) {
    .Deprecated("searchProtein")
    return searchProtein(pattern)
}


getAnnotation <- function(query){
    .Deprecated("annotateSequence")
    return annotateSequence(query)
}

getGenomeAlignment <- function(genome1, genome2, chr1, chr2){
    .Deprecated("getGenomePair")
    return getGenomePair(genome1, genome2, chr1, chr2)
}

getOntologies <- function(df){
    .Deprecated("getProtein")
    return(getProtein(ids=df$oma_id)$gene_ontology)
}

getSequences <- function(df){
    .Deprecated("getProtein")
    return(getProtein(ids=df$oma_id)$sequence)
}

getData <- function(type, id, attribute){
    if (type == "group"){
        .Deprecated("getOMAGroup")
        return(getOMAGroup(id, attribute))
    } else if (type == "genome"){
        .Deprecated("getGenome")
        return(getGenome(id, attribute))
    } else if (type == "protein"){
        .Deprecated("getProtein")
        return(getProtein(id, attribute))
    } else {
        .Deprecated("getProtein")
        return(NULL)
    }
}
