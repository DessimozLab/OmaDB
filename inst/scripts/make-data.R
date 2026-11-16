## Script to generate the example datasets in data/.
##
## All datasets are responses of the OMA Browser REST API
## (https://omabrowser.org/api/) retrieved with the functions of this
## package. Run this script from the root of the package source with the
## package installed (or loaded with devtools::load_all()). The data
## corresponds to the OMA release reported by getVersion() at the time of
## generation (All.May2026).

library(OmaDB)

seq <- paste0("MNDPSLLGYPNVGPQQQQQQQQQQHAGLLGKGTPNALQQQLHMNQLTGIPPPGLMN",
    "NSDVHTSSNNNSRQLLDQLANGNANMLNMNMDNNNNNNNNNNNNNNNGGGSGVMMNASTAAVNSIGM")

protein <- getProtein("YEAST00001")
orthologs <- getProtein("YEAST00001", attribute = "orthologs")
group <- getOMAGroup("YEAST00058")
hog <- getHOG("YEAST00001")
pairs <- utils::head(getGenomePairs("YEAST", "ASHGO", chr1 = "I"), 100)
taxonomy <- getTaxonomy(members = c("YEAST", "ASHGO", "SCHPO"))
xref <- searchProtein("P53_RAT")
sequence_map <- mapSequence(seq)
sequence_annotation <- annotateSequence(seq)

for (name in c("protein", "orthologs", "group", "hog", "pairs",
        "taxonomy", "xref", "sequence_map", "sequence_annotation")) {
    save(list = name, file = file.path("data", paste0(name, ".rda")),
        compress = "xz")
}
