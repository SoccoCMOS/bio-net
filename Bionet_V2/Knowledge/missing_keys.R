all_names<-read.csv2("0.full_gbif_normalized_names.csv")
all_names$occurrenceId<-1:nrow(all_names)
missing_names<-read.csv2("taxons_compl.csv")
missing_names$occurrenceId<-nrow(all_names)+1:nrow(missing_names)
full<-unique(rbind(all_names,missing_names))
still_missing<-subset(full,matchType=="NONE")$verbatimScientificName

#### Get keys from troph
troph<-read.csv2("1.full_interactions_list.csv")
nm_cons=troph[c("consumer_gbif_key","consumer_name")]
colnames(nm_cons)=c("key","verbatimScientificName")
nm_res=troph[c("resource_gbif_key","resource_name")]
colnames(nm_res)=c("key","verbatimScientificName")
nm_troph=unique(rbind(nm_cons,nm_res))


#### Check conflicts inside gbif
#### Select most precise assignation for taxa with different assignations
conflict_gbif=list()
### Remove NA keys (those with no assignation in GBIF)
full_notna=subset(full,!is.na(key))
full_notna$key=as.character(full_notna$key)

### Treat those with multiple assignations
all_freq=as.data.frame(table(full_notna$key))
for (name in full_notna$verbatimScientificName){
  sub=subset(full_notna,verbatimScientificName==name)
  if(length(unique(sub$key))>1){
    conflict_gbif=append(conflict_gbif,name)
    freq=all_freq[all_freq$Var1%in%sub$key,]
    new_key=freq[which.min(freq$Freq),]$Var1
    
    ####Find a way to drop the instances that do not match this key
    full_notna=subset(full_notna,(verbatimScientificName==name & key==new_key) | (verbatimScientificName!=name))
  }
}

###Keep the ones that have been matched when conflict
conflict_gbif_post=list()
all_freq_post=as.data.frame(table(full_notna$key))
for (name in full_notna$verbatimScientificName){
  sub=subset(full_notna,verbatimScientificName==name)
  if(nrow(sub)>1){
    conflict_gbif_post=append(conflict_gbif_post,name)
    full_notna=subset(full_notna,(verbatimScientificName==name & !is.na(scientificName) & !is.null(scientificName)) | (verbatimScientificName!=name))
  }
}

###Drop duplicates verbatim with same info
clean_assign=unique(full_notna[,-1])

###TODO the problem of the duplicate keys => apply post-coding procedure
###TODO the problem of the NA keys => drop 
###TODO: correct keys in troph from clean_assign


write.csv2(full,"0.full_taxa.csv",row.names = FALSE)

