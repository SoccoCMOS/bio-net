
proba_int<-read.csv2("50__groups.csv")[,-1]
row.names(proba_int)<-colnames(proba_int)

label_assign <- read.csv2("50__labels.csv")
taxonomy <- read.csv2("../Knowledge/0.taxonomy_augmented.csv")

##Augment label assign with rank and name
cluster_augm<-merge(label_assign,taxonomy,by.x = "X",by.y = "key","join")[c("X","verbatimScientificName","rank","labels","kingdom","phylum","class","order","family")]

write.csv2(cluster_augm,"clusters.csv",row.names = F)
