library(dplyr)
###Check taxo and names and keys
gbif_names <- read.csv2("0.gbif_normalized_names.csv",encoding = "latin-1")

###OccurrenceIDs
gbif_names$occurrenceId=1:dim(gbif_names)[1]

###Duplicate keys
row.names(gbif_names)<-gbif_names$key

duplicates_returned<-c("48", "2120985", "2121983", "2193785", "2282644", "2284116", "2308266", "2308278", "2308378", "2638548", "2639987", "3194385", "3194398", "3194566", "3200900", "3201016", "3201179", "3201190", "3201319", "3202401", "3210222", "3224243", "3236919", "3285828", "4281960", "4378440", "4378507", "4892213", "4893086", "4899951", "4900166", "5271044", "5899695", "7073608", "7566978", "7575131", "7757178", "7787742", "7990941", "8045705", "8079935", "8167254", "8198844", "8280877", "8410027", "9171551", "9381748", "9530095", "9762970")

###Returns synonyms <- same taxo different verbatimScientificName 
gbif_dupl<-subset(gbif_names,key%in%duplicates_returned)
syn=vector(mode="list",length(duplicates_returned))
cpt=0
for (k in duplicates_returned){
  cpt=cpt+1
  syn[[cpt]]$k=k
  syn[[cpt]]$syn=as.character(subset(gbif_dupl,key==k)$verbatimScientificName)
}

###Unknown keys
ukn=as.character(subset(gbif_names,is.na(key))$verbatimScientificName)
write.csv2(ukn,"unknown_keys.csv")

###Check interaction list
int_list <- read.csv2("1.interactions_list.csv",encoding = "latin-1")

##get keys and names for consumer+resource
cons=int_list[c("consumer_gbif_key","consumer_name")]
colnames(cons)=c("key","name")
res=int_list[c("resource_gbif_key","resource_name")]
colnames(res)=c("key","name")

all_names=distinct(rbind(cons,res))
rownames(all_names)<-all_names$key

duplicate_keys=c("1", "7", "2410", "3565", "2704913", "2873815", "2927173", "2973363", "7995850", "8395064" )
check_duplicates=subset(all_names,key%in%duplicate_keys)

###Merge with verified
merged<-subset(all_names,key %in% gbif_names$key)