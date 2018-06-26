setwd("J:/PROJET RECONSTRUCTION DE RESEAUX/TEST MATRICE TUTEUR/Script 21062018/Comparaison clustering")
install.packages('clusteval', dependencies = TRUE)
library(clusteval)

#Importation des data : taxa + classes d'appartenance
SBM          <- read.csv("SBM_Poiss_class.csv", sep=";", header = T)
hierarchical <- read.csv("groupsW-10.csv", sep=";", header=T)
spectral     <- read.csv("10_communities.csv", sep=";", header=T)

#Ajout des familles et du code FW
metabar     <- read.csv("mb.csv", sep = ";", h = T) 

label_SBM  <- sort(SBM_clust$retained_tax) #tri les noms de taxa dans l'ordre alphabétique
label_hier <- sort(hierarchical_clust$retained_tax)
label_spec <- sort(spectral_clust$SPEC)

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

#Tableau complet : taxa, class, code FW, fam
SBM_clust          <- read.csv("label_SBM.csv", sep=";", header = T)
hierarchical_clust <- read.csv("label_spec.csv", sep=";", header=T)
spectral_clust     <- read.csv("label_hier.csv", sep=";", header=T)

# Comparaison entre le modèle SBM et la classification experte
eval <- comembership_table(SBM_clust$class,SBM_clust$codeFW)

# Comparaison entre le modèle hiérarchique et la classification experte
eval2 <- comembership_table(hierarchical_clust$class,hierarchical_clust$codeFW)

# Comparaison entre le modèle spectral et la classification experte
eval3 <- comembership_table(spectral_clust$class,spectral_clust$codeFW)

#Comparaison des modèles deux à deux
hier_spec <- comembership_table(hierarchical_clust$class, spectral_clust$LABEL)
SBM_spec  <- comembership_table(SBM_clust$class, spectral_clust$LABEL)
SBM_hier  <- comembership_table(SBM_clust$class, hierarchical_clust$class)

#Similarité entre les modèles comparés
#Calcul suivant un indice de Jaccard tel que J = 11/(11+10+01)
similar_hier_spec <- hier_spec$n_11/(hier_spec$n_11 + hier_spec$n_10 + hier_spec$n_01)
similar_SBM_spec <- SBM_spec$n_11/(SBM_spec$n_11 + SBM_spec$n_10 + SBM_spec$n_01)
similar_SBM_hier <- SBM_hier$n_11/(SBM_hier$n_11 + SBM_hier$n_10 + SBM_hier$n_01)

#Verif possible
#cluster_similarity(hierarchical_clust$class, spectral_clust$LABEL, similarity = c("jaccard"), method = "independence")




