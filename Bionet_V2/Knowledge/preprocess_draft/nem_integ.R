nem_names<-read.csv("nematodes/gbif_species_names.csv")[,1:11]
all_names<-read.csv2("0.gbif_normalized_names.csv")
  
nem_troph<-read.csv("nematodes/trophic_group.csv")
all_troph<-read.csv2("1.interactions_list.csv")

nem_size<-read.csv("nematodes/body_length.csv")
all_size<-read.csv2("2.body_length.csv")

###1. Remove old info on trophic groups
drop_nem<-all_troph[!(all_troph$consumer_gbif_key %in% nem_troph$consumer_gbif_key),]

####2. Add new info
##Align data structures
cols<-names(all_troph)
insert_nem<-nem_troph[cols]
new_troph<-rbind(drop_nem,insert_nem)

write.csv2(new_troph,"1.full_interactions_list.csv",row.names = FALSE)

###3. Remove old info on body size
drop_nemsize<-all_size[!(all_size$gbif_key %in% nem_size$gbif_key),]

###4. Add new info
colsize<-names(all_size)
insert_nemsize<-nem_size[colsize]
new_size<-rbind(drop_nemsize,insert_nemsize)
write.csv2(new_size,"2.full_body_length.csv",row.names = FALSE)