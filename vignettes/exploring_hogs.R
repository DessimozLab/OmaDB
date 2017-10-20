## ---- warning=FALSE, message=FALSE---------------------------------------
library(roma)

hog <- getHOG(id="HUMAN22168")

getObjectAttributes(hog)

hog_id = getAttribute(hog,'hog_id')

parent_hogs = getAttribute(hog,'parent_hogs')

dim(parent_hogs)

parent_hog_id = parent_hogs[[hog_id]]

children_hogs = getAttribute(hog,'children_hogs')

children_hogs

child_hog_ids  = children_hogs[['hog_id']]

child_hog_1 = getHOG(id=child_hog_ids[[1]])

getAttribute(child_hog_1,"level")



## ---- warning=FALSE, message=FALSE---------------------------------------

childrenMembers = childrenMembers(hog)

myInterestingGenes = childrenMembers[['HOG:0261495.1a.1a']]

annotations = getOntologies(childrenMembers)

# We can now directly use the above data to construct a topGO object for further analysis (such as GO enrichment)

library(topGO)

topGO = getTopGO(annotations, myInterestingGenes, format = "geneID2GO")



