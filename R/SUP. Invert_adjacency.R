setwd("C:\Users\simoussi\ownCloud\Ma_these\bio-net\data\input\taxo_metabar\adj\sbm")
spec=read.csv2("species.csv",row.names=1,header=T)
family=read.csv2("family.csv",row.names=1,header=T)
genus=read.csv2("genus.csv",row.names=1,header=T)


spec_t=data.frame(t(spec))

family_t=data.frame(t(family))

genus_t=data.frame(t(genus))

write.csv2(spec_t,"t/species.csv")
write.csv2(genus_t,"t/genus.csv")
write.csv2(family_t,"t/family.csv")
