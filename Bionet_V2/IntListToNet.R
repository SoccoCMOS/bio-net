library(dplyr)
library(rlist)
library(reshape2)

int_list<-read.csv2("Knowledge/TrophParasInteractions.csv")
taxa_list <- read.csv2("Data/normalized_metabar_out.csv")
taxonomy <- read.csv2("Knowledge/0.taxonomy_augmented.csv")

###Removing empty spaces from taxa names on metabar file
taxa_list$verbatimScientificName<-as.character(lapply(taxa_list$verbatimScientificName,function(x) strTrim(as.character(x))))

###Merge info from taxonomy for each MOTU: key, rank
taxa_list_compl<-merge(taxa_list,taxonomy[c("key","verbatimScientificName","rank")],by = "verbatimScientificName")

###Subselect only interactions involving taxa from our pool
sub_int=subset(int_list,consumer %in% taxa_list_compl$key|resource %in% taxa_list_compl$key )
write.csv2(sub_int,"Adjacency_list.csv",row.names = FALSE)

###Create network