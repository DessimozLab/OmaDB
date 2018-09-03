## ---- warning=FALSE, message=FALSE---------------------------------------
library(OmaDB)

xref = load('../data/xref.rda')

head(xref)

## ---- warning=FALSE, message=FALSE---------------------------------------

load('../data/pairs.rda')

head(pairs)


## ---- warning=FALSE, message=FALSE---------------------------------------

load('../data/group.rda')

object_attributes = getObjectAttributes(group)

group$fingerprint

getAttribute(group, 'fingerprint')



## ---- warning=FALSE, message=FALSE---------------------------------------

load('../data/protein.rda')

getAttribute(protein,'orthologs')

load('../data/orthologs.rda')

orthologs 


## ---- warning=FALSE, message=FALSE---------------------------------------

gRanges = getInfo(orthologs,type='genomic_ranges')

str(gRanges)


