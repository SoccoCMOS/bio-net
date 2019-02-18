####Read data 
normalized_names=read.csv2("1.normalized.csv")
mb_list=read.csv2("metabar_out.csv")

####Normalize MOTUs assignations
mb_list$verbatimScientificName=as.character(mb_list$scientific_name)
for(i in 1:nrow(normalized_names)){
  old=as.character(normalized_names[i,"verbatimScientificName"])
  norm=as.character(normalized_names[i,"retainedName"])
  nb=nrow(mb_list[mb_list$verbatimScientificName==old & !is.na(mb_list$verbatimScientificName),])
  if(nb>0)
    mb_list[mb_list$verbatimScientificName==old & !is.na(mb_list$verbatimScientificName),"verbatimScientificName"]<-rep(norm,nb)
}

write.csv2(mb_list,"normalized_mb.csv",row.names=FALSE)
