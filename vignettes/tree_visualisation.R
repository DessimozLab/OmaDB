## ---- warning=FALSE, message=FALSE---------------------------------------
library(roma)

taxonomy = getTaxonomy(root="Alveolata")

write.table(taxonomy$newick, "newick.txt")

tree_object = getTree("newick.txt")

ggtree::ggtree(tree_object)



