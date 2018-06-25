# Working directory

el          <- read.csv("input/expert_edge_list.csv", sep = ";", h = T) #TG1, TG2 et type d'interaction : symbiotique, trophique, ou parasitique
metabar     <- read.csv("input/metabar.csv", sep = ";", h = T) 

TG1 <- el$TG1
TG2 <- el$TG2
type <- union(TG1,TG2)
codeFW <- metabar$codeFW
setdiff(type, codeFW) 
setdiff(codeFW, type)
