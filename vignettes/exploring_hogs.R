## ---- warning=FALSE, message=FALSE---------------------------------------
library(roma)

hog <- getHOG(id="HUMAN22168")

getObjectAttributes(hog)

hog_id = getAttribute(hog,'hog_id')

parent_hogs = getAttribute(hog,'parent_hogs')

dim(parent_hogs)

parent_hog_id = parent_hogs$hog_id[[1]]

children_hogs = getAttribute(hog,'children_hogs')

child_hog_id_1  = children_hogs$hog_id[[1]]
child_hog_id_1

child_hog_id_2  = children_hogs$hog_id[[2]]
child_hog_id_2

child_hog_1 = getHOG(id=child_hog_id_1)

getAttribute(child_hog_1,"level")



## ---- warning=FALSE, message=FALSE---------------------------------------
library(topGO)

child_hog_1_members = getHOG(id = child_hog_id_1, members=TRUE)
child_hog_2_members = getHOG(id = child_hog_id_2, members=TRUE)

# We can now directly use the above data to construct a topGO object for further analysis (such as GO enrichment)

annotations_1 = getOntologies(child_hog_1_members$members)
annotations_2 = getOntologies(child_hog_2_members$members)

annotations = c(annotations_1,annotations_2)

topGO = getTopGO(annotations, myInterestingGenes = names(annotations_1), format = "geneID2GO")



