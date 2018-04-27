## ---- warning=FALSE, message=FALSE---------------------------------------
library(roma)

load('../data/sequence_map.rda')

getObjectAttributes(sequence_map)

targets = getAttribute(sequence_map,'targets')

length(targets) 

protein = targets[[1]][['entry_url']]


## ---- warning=FALSE, message=FALSE---------------------------------------

load('../data/sequence_annotation.rda')

sequence_annotation 


