
names<-read.csv2("0.Completed_taxa_list.csv")
interactions<-read.csv2("1.Completed_interaction_list.csv",stringsAsFactors = FALSE)
body_length<-read.csv2("2.full_body_length.csv")
micro_hab<-read.csv2("3.micro_hab.csv")

###Complete interactions keys from names keys
no_cons_key=subset(interactions,is.na(consumer_gbif_key))$consumer_name
no_res_key=subset(interactions,is.na(resource_gbif_key))$resource_name

no_key=union(no_cons_key,no_res_key)

get_key<-function(n){
  n=trimws(as.character(n))
  return(as.character(subset(names,verbatimScientificName==n | scientificName==n)$key))
}

keys=lapply(no_key, get_key)

df_miss=data.frame(no_key)
df_miss$key=keys

####Refill missing keys
for (m in no_cons_key){
  interactions[interactions$consumer_name==m,"consumer_gbif_key"]=get_key(m)
}

for (m in no_res_key){
  interactions[interactions$resource_name==m,"resource_gbif_key"]=get_key(m)
}

write.csv2(interactions,"1.interactions_list.csv")

names=read.csv2("0.Completed_taxa_list.csv")
int=read.csv2("1.interactions_list.csv")

diff_res=setdiff(int$resource_name,names$verbatimScientificName)
diff_cons=setdiff(int$consumer_name,names$verbatimScientificName)

int$consumer_name=lapply(int$consumer_name,trimws)
int$resource_name=lapply(int$resource_name,trimws)