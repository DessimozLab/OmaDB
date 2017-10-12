## ---- warning=FALSE, message=FALSE---------------------------------------
library(roma)

hog <- getHOG(id="HUMAN22168")

attributes(hog)

hog$hog_id

length(hog$parent_hogs)

hog$parent_hogs[[1]]$parent_hogs

length(hog$children_hogs)

hog$children_hogs[[1]]$hog_id

hog$children_hogs[[2]]$hog_id

roma::resolveURL(hog$children_hogs[[1]]$levels_url)$level



## ---- warning=FALSE, message=FALSE---------------------------------------
library(roma)
library(topGO)

hog <- getHOG(id="HUMAN22168")

children = hog$children_hogs


#creating a list of omaid's for the protein members
proteins_list = list()

for(child in children){
	child_hog = resolveURL(child$levels_url)
	members = resolveURL(child_hog$members_url)$members
	for(i in range(length(members))){
		proteins_list[[i]] = members[[i]]$omaid
	}
}

#let's check the go annotations for each

object_list = list()

for(protein in proteins_list){
	object_list[[protein]] = getData(type="protein",id=protein)
}

# We can now directly use the above data to construct a topGO object for further analysis (such as GO enrichment)

annotations = formatTopGO(object_list,format="geneID2GO")

topGO = getTopGO(annotations, myInterestingGenes = list(proteins_list[[1]]), format = "geneID2GO")



