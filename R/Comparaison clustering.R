setwd("J:/PROJET RECONSTRUCTION DE RESEAUX/TEST MATRICE TUTEUR/Script 21062018/Comparaison clustering")
install.packages('clusteval', dependencies = TRUE)
library(clusteval)

SBM_clust <- read.csv("SBM_Poiss_class.csv", sep=";", header = T)

eval<-comembership_table(SBM_clust$class,SBM_clust$codeFW)

hierarchical_clust <- read.csv("groupsW-10.csv", sep=";", header=T)
spectral_clust <- read.csv("10_communities.csv", sep=";", header=T)

label_SBM  <- sort(SBM_clust$retained_tax)
label_hier <- sort(hierarchical_clust$retained_tax)
label_spec <- sort(spectral_clust$SPEC)
t <- comembership_table(hierarchical_clust$class, spectral_clust$LABEL)
u <- comembership_table(SBM_clust$class, spectral_clust$LABEL)
v <- comembership_table(SBM_clust$class, hierarchical_clust$class)

#Ajout des familles et du code FW
metabar     <- read.csv("mb.csv", sep = ";", h = T) 

label_hier$fam <- NULL
label_hier$codeFW <- NULL
for (i in 1:nrow(label_hier)) {
  sp <- label_hier$retained_tax[i]
  sp <- as.character(sp[1])
  fam <- unique(metabar$family_name[which(metabar$retained_tax == sp)])
  label_hier$fam[i] <- as.character(fam[1])
  FW <- unique(metabar$codeFW[which(metabar$retained_tax == sp)])
  label_hier$codeFW[i] <- as.character(FW[1])
}
write.csv2(label_hier, "label_hier.csv")

label_SBM$fam <- NULL
label_SBM$codeFW <- NULL
for (i in 1:nrow(label_SBM)) {
  sp <- label_SBM$retained_tax[i]
  sp <- as.character(sp[1])
  fam <- unique(metabar$family_name[which(metabar$retained_tax == sp)])
  label_SBM$fam[i] <- as.character(fam[1])
  FW <- unique(metabar$codeFW[which(metabar$retained_tax == sp)])
  label_SBM$codeFW[i] <- as.character(FW[1])
}
write.csv2(label_SBM, "label_SBM.csv")

label_spec$fam <- NULL
label_spec$codeFW <- NULL
for (i in 1:nrow(label_spec)) {
  sp <- label_spec$SPEC[i]
  sp <- as.character(sp[1])
  fam <- unique(metabar$family_name[which(metabar$retained_tax == sp)])
  label_spec$fam[i] <- as.character(fam[1])
  FW <- unique(metabar$codeFW[which(metabar$retained_tax == sp)])
  label_spec$codeFW[i] <- as.character(FW[1])
}
write.csv2(label_spec, "label_spec.csv")