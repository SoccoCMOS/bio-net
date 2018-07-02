# -----------------------------------------------------------------------------
# PROJECT:
#    BISE : Linking soil foodweb and soil functions dynamics under
#    agricultural practices (des)intensification and different input levels
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# STEP:
#    0.1   Generating species x species adjacency matrices (binary, weighted by proba & count)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# NOTES:
#  
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# LIBRARIES:
library(reshape2)
# -----------------------------------------------------------------------------

##### Params ####

filt_type_inter=c("symbiotic")
#################

setwd("../data/input/taxo_metabar")

# -----------------------------------------------------------------------------
# data load
# -----------------------------------------------------------------------------
#Metabarcoding dataset
metabar     <- read.csv("family_family_metabar_raw.csv", sep = ";", h = T, row.names=1) 
dim(metabar)
#edge list of expert based interactions
el          <- read.csv("../expert_edge_list.csv", sep = ";", h = T)   ##Expert edges
el         <- el[!(el$type %in% filt_type_inter),]  ### Filter specific types of interactions
dim(el)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Modifications on metabarcoding dataframe
# -----------------------------------------------------------------------------
# getting a dataframe with only taxa name, their expert-based trophic group and their µhabitat => expert_memb
mb <- metabar[,c("retained_tax", "codeFW", "µhab_surf", "µhab_subsurf", "µhab_soil")]
expert_memb <- mb[complete.cases(mb),] ##Expert groupings, remove NAN
expert_memb <- unique(expert_memb)  ##Keep unique values

full_taxa <- unique(expert_memb$retained_tax) ### Full list of retained taxa
length(full_taxa)

polyv_taxa <- expert_memb$retained_tax[which(duplicated(expert_memb$retained_tax))]
print(list(polyv_taxa))
length(polyv_taxa)

# getting a dataframe with taxa presence/absence by sampling point  => metabar_bin
metabar_agg <- aggregate(metabar[, -c(1:16,273)], list(retained_tax = metabar$retained_tax), sum, na.rm=TRUE) # sum by taxon name
rownames(metabar_agg) <- metabar_agg$retained_tax
metabar_agg <- metabar_agg[, -1]
metabar_bin <- ifelse(metabar_agg == 0, 0, 1) ### => Presence/Absence of retained taxa per plot (BISExy)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Building a binary species adjacency matrix 
# -----------------------------------------------------------------------------
# building a "species x species" edge list from a "trophic group x trophic group" edge list
el2 <- merge(x = el, y = expert_memb, by.x = "TG2", by.y = "codeFW", all = FALSE)
el3 <- merge(x = el2, y = expert_memb, by.x = "TG1", by.y = "codeFW", suffixes = c("_TG2", "_TG1"))

#### Number of species after merge = Species that have targeted (non fitlered) interactions
interac_taxa <- union(el3$retained_tax_TG1,el3$retained_tax_TG2)
unknown<-setdiff(full_taxa,interac_taxa) ### Species with no background knowledge on interaction behavior or species only involved in filtered interaction types as of expert knowledge
length(interac_taxa)
length(unknown)
####

#fitering by µhabitat
el3$prof_match   <- el3$µhab_surf_TG1 * el3$µhab_surf_TG2 + ##Co-occurrence at the surface level
  el3$µhab_subsurf_TG1 * el3$µhab_subsurf_TG2+ ### OR co-occurrence at the subsurface level
  el3$µhab_soil_TG1 * el3$µhab_soil_TG2   ### OR co-occurrence at the soil level

el3$interac   <- ifelse(el3$prof_match == 0, 0, 1) ### Interac is true if co-occurrence match in at least one µhabitat, false else.

spec_adj      <- dcast(data = el3, retained_tax_TG1~retained_tax_TG2, value.var = "interac", fill = 0, fun.aggregate = sum, drop = F) ### spec_adj_bin= adjacency matrix for taxa-taxa interactions
rownames(spec_adj)<- spec_adj[,1] ## Index rows by taxa name
spec_adj_sq <- spec_adj[,-1] ## Remove first column => obtain a square adjacency matrix

# deleting lines and columns for which there is no interactions
condition <- as.data.frame(rowSums(spec_adj_sq))
condition$colSums <- colSums(spec_adj_sq)
condition$cond <- rowSums(condition)
spec_adj_bin <- as.matrix(spec_adj_sq[condition$cond  > 0, condition$cond > 0])
dim(spec_adj_bin)
others<-setdiff(interac_taxa,rownames(spec_adj_bin))  ### Other taxa filtered by µhabitat
others
length(others)

# -----------------------------------------------------------------------------
# Building a weighted species adjacency matrix 
# -----------------------------------------------------------------------------
#assessing the probability of interaction under a neutral hypothesis
mbf              <- as.matrix(subset(metabar_bin, rownames(metabar_bin) %in% rownames(spec_adj_bin)))
cooccur              <- mbf %*% t(mbf) ### Co-occurrence counts within plots between every pair of species
out              <- cooccur/256 ### Probability of co-occurrence = counts of plots of cooccurrence / number of plots in total

dim(out)

spec_adj_weighted         <- out*spec_adj_bin  ### Weighted interaction matrix (following neutral hypothesis,  weight=co-occurrence probability)
spec_adj_count            <- cooccur*spec_adj_bin  ### Replace co-occurrence probability by count
dim(spec_adj_weighted)
dim(spec_adj_count)
dim(cooccur)

#### Save adjacency matrices ###
# write.csv2(metabar_bin,"metabar_bin.csv")
# write.csv2(mbf,"filtered_metabar.csv")
# write.csv2(t(mbf),"transpose_filtered_metabar.csv")
# write.csv2(cooccur,"cooccurrence_proba.csv")
write.csv2(spec_adj_bin, "adjacence_esp_bin_family.csv")
# write.csv2(spec_adj_weighted, "adjacence_esp_weighted.csv")
# write.csv2(spec_adj_count, "adjacence_esp_count.csv")
# -----------------------------------------------------------------------------
