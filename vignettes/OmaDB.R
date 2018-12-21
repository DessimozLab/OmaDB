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



