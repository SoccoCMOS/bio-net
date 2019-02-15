all_names<-read.csv2("0.full_gbif_normalized_names.csv")
all_names$occurrenceId<-1:nrow(all_names)
completed_names<-read.csv2("taxons_compl.csv")
completed_names$occurrenceId<-nrow(all_names)+1:nrow(completed_names)
full<-unique(rbind(all_names,completed_names))

still_missing<-unique(data.frame(ScientificName=subset(full,matchType=="NONE" | is.na(key))$verbatimScientificName))

#### Get keys from troph
troph<-read.csv2("1.full_interactions_list.csv")

#### Gather all names and keys to complete them from GBIF
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
clean_assign=unique(full_notna)
all_freq_clean=as.data.frame(table(clean_assign$key))

###Process keys from troph
##Find unassigned taxa
nm_troph_unass=unique(subset(nm_troph,!(verbatimScientificName %in% all_names$verbatimScientificName)))
write.csv2(nm_troph_unass,"troph_unassigned.csv",row.names = FALSE)
##Read assignation (result of GBIF query) for taxa that are not covered in assign => 1130-1113=17 are not found anywhere
troph_compl_assign<-read.csv2("troph_unass_compl.csv")
full_assign<-unique(rbind(clean_assign,troph_compl_assign))

###GOTO check conflicts procedure
full_notna=full_assign

write.csv2(clean_assign,"0.Completed_taxa_list.csv",row.names = FALSE)

full_assign=clean_assign

###Correct wrong names in troph then merge keys to troph
correct_names<-read.csv2("correct_names.csv")
correct_tax_name<-function(n){
  n=as.character(n)
  if(n%in%correct_names$oldScientificName)
    return(as.character(correct_names[correct_names$oldScientificName==n,"verbatimScientificName"]))
  else
    return(n)
}

get_key<-function(n){
  sub=unique(subset(full_assign,verbatimScientificName==n))
  keys=sub$key
  ranks=sub$rank
  
  if(length(keys)>1)
    print(n)
  return(list("k"=as.character(keys[1]),"r"=as.character(ranks[1])))
}

treat_troph<-function(x){
  ###Correct name eventually
  x$resource_name<-correct_tax_name(x$resource_name)
  x$consumer_name<-correct_tax_name(x$consumer_name)
  ###Check key and rank 
  res=get_key(x$resource_name)
  cons=get_key(x$consumer_name)
  
  x$resource_gbif_key<-res$k 
  x$resource_rank<-res$r
  
  x$consumer_gbif_key<-cons$k
  x$consumer_rank<-cons$r
  
  return(x)
}

#troph_corrected<-apply(troph,1,treat_troph)
new_troph=data.frame(matrix(ncol = ncol(troph), nrow = 0))
colnames(new_troph) <- colnames(troph)

for (i in 1:nrow(troph)){
  x<-troph[i,]
  new_x<-treat_troph(x)
  new_troph=rbind(new_troph,new_x)
}
#interactions<-apply(troph_corrected,2,as.character)

###TODO the problem of the duplicate keys on full_assign => apply post-coding procedure
write.csv2(new_troph,"1.Completed_interaction_list.csv",row.names = FALSE)

still_missing=rbind(still_missing,data.frame("ScientificName"=subset(new_troph,is.na(consumer_gbif_key))$consumer_name))
still_missing=rbind(still_missing,data.frame("ScientificName"=subset(new_troph,is.na(resource_gbif_key))$resource_name))
write.csv2(still_missing,"not_found.csv",row.names = FALSE)
