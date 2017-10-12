## ---- warning=FALSE, message=FALSE---------------------------------------
library(roma)

results <- getXref(pattern="MAL")

results

## ---- warning=FALSE, message=FALSE---------------------------------------
library(roma)

alignment <- getGenomeAlignment(genome_id1 = "YEAST", genome_id2 = "ASHGO")

head(alignment)


## ---- warning=FALSE, message=FALSE---------------------------------------
library(roma)

group <- getData(type="Group",id="YEAST58")

group$fingerprint



## ---- warning=FALSE, message=FALSE---------------------------------------
library(roma)

protein <- getData(type="Protein",id="YEAST58")

attributes(protein)

protein$orthologs

orthologs = resolveURL(protein$orthologs)

head(orthologs)

formatted_orthologs = formatData(orthologs)


