## ---- warning=FALSE, message=FALSE---------------------------------------
library(roma)

results <- getXref(pattern="MAL")

head(results)

## ---- warning=FALSE, message=FALSE---------------------------------------

alignment <- getGenomeAlignment(genome_id1 = "YEAST", genome_id2 = "ASHGO")

head(alignment)


## ---- warning=FALSE, message=FALSE---------------------------------------

group <- getData(type="group",id="YEAST58")

group$fingerprint

getAttribute(group, 'fingerprint')



## ---- warning=FALSE, message=FALSE---------------------------------------

protein <- getData(type="protein",id="YEAST58")

attributes(protein)

protein$orthologs

orthologs = resolveURL(protein$orthologs)


