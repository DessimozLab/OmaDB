## ---- warning=FALSE, message=FALSE---------------------------------------
library(OmaDB)

load('../data/hog.rda')

print(hog)

getObjectAttributes(hog)

hog_id = getAttribute(hog,'hog_id')

parent_hogs = getAttribute(hog,'parent_hogs')

dim(parent_hogs)

parent_hog_id = parent_hogs[[hog_id]]

children_hogs = getAttribute(hog,'children_hogs')

children_hogs

child_hog_ids  = children_hogs[['hog_id']]



