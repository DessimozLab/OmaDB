## ---- warning=FALSE, message=FALSE---------------------------------------
library(roma)

taxonomy = getTaxonomy(root="Alveolata")

tree_object = getTree(taxonomy)

ggtree::ggtree(tree_object)



