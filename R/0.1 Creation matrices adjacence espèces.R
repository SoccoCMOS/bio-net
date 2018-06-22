# -----------------------------------------------------------------------------
# PROJECT:
#    BISE : Linking soil foodweb and soil functions dynamics under
#    agricultural practices (des)intensification and different input levels
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# STEP:
#    0.1   Generating species x species adjacency matrices
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# NOTES:
#  
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# LIBRARIES:
library(reshape2)
# -----------------------------------------------------------------------------


setwd("C:/1. Mike/5. RH/1. Stagiaire/18_Labouyrie Maeva/v2/bio-net/Data")

# -----------------------------------------------------------------------------
# data load
# -----------------------------------------------------------------------------
#Metabarcoding dataset
metabar     <- read.csv("mb.csv", sep = ";", h = T) 
#edge list of expert based interactions
el          <- read.csv("expert_edge_list.csv", sep = ";", h = T)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Modifications on metabarcoding dataframe
# -----------------------------------------------------------------------------
# getting a dataframe with only taxa name, their expert-based trophic group and their µhabitat => expert_memb
mb <- metabar[,c("retained_tax", "codeFW", "µhab_surf", "µhab_subsurf", "µhab_soil")]
expert_memb <- mb[complete.cases(mb),]
expert_memb <- unique(expert_memb) 
# getting a dataframe with taxa presence/absence by sampling point  => metabar_bin
metabar_agg <- aggregate(metabar[, -c(1:17)], list(retained_tax = metabar$retained_tax), sum) # sum by taxon name
rownames(metabar_agg) <- metabar_agg$retained_tax
metabar_agg <- metabar_agg[, -1]
metabar_bin <- ifelse(metabar_agg == 0, 0, 1)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Building a binary species adjacency matrix 
# -----------------------------------------------------------------------------
# building a "species x species" edge list from a "trophic group x trophic group" edge list
el2 <- merge(x = el, y = expert_memb, by.x = "TG2", by.y = "codeFW", all = FALSE)
el3 <- merge(x = el2, y = expert_memb, by.x = "TG1", by.y = "codeFW", suffixes = c("_TG2", "_TG1"))

#fitering by µhabitat
el3$prof_match   <- el3$µhab_surf_TG1 * el3$µhab_surf_TG2 + 
  el3$µhab_subsurf_TG1 * el3$µhab_subsurf_TG2+
  el3$µhab_soil_TG1 * el3$µhab_soil_TG2
el3$prof_proba   <- ifelse(el3$prof_match == 0, 0, 1)
el4               <- el3[el3$type != "symbiotic",]
spec_adj_bin      <- dcast(data = el4, retained_tax_TG1~retained_tax_TG2, value.var = "prof_proba", fill = 0, fun.aggregate = sum, drop = F)
rownames(spec_adj_bin)<- spec_adj_bin[,1]
spec_adj_bin <- spec_adj_bin[,-1]

# deleting lines and columns for which there is no interactions
condition <- as.data.frame(rowSums(spec_adj_bin))
condition$colSums <- colSums(spec_adj_bin)
condition$cond <- rowSums(condition)
spec_adj_bin <- spec_adj_bin[condition$cond  > 0, condition$cond > 0]
spec_adj_bin <- as.matrix(spec_adj_bin)
write.csv2(spec_adj_bin, "adjacence_esp_bin.csv")

# -----------------------------------------------------------------------------
# Building a weighted species adjacency matrix 
# -----------------------------------------------------------------------------
#assessing the probability of interaction under a neutral hypothesis
num              <- (metabar_bin %*% t(metabar_bin))
out              <- num/256
neutral          <- melt(out, value.name = "neutral_prob")
neutral$code     <- paste(neutral$Var1, neutral$Var2, sep = "_")
el3$code         <- paste(el3$retained_tax_TG1, el3$retained_tax_TG2, sep = "_")
el3              <- merge(el3, neutral, "code")
el3$proba        <- el3$prof_proba * el3$neutral_prob
el3              <- el3[!is.na(el3$proba),]

# bilding the species adjacency matrix without symbiotic interactions
el3                        <- el3[el3$type != "symbiotic",]
spec_adj_weighted          <- dcast(data = el3, retained_tax_TG1~retained_tax_TG2, value.var = "proba", fun.aggregate = sum, drop = FALSE)
rownames(spec_adj_weighted)<- spec_adj_weighted[,1]
spec_adj_weighted          <- spec_adj_weighted[,-1]

# deleting lines and columns for which there is no interactions
condition <- as.data.frame(rowSums(spec_adj_weighted))
condition$colSums <- colSums(spec_adj_weighted)
condition$cond <- rowSums(condition)
spec_adj_weighted <- spec_adj_weighted[condition$cond  > 0, condition$cond > 0]
spec_adj_weighted <- as.matrix(spec_adj_weighted)
write.csv2(spec_adj_weighted, "adjacence_esp_weighted.csv")
# -----------------------------------------------------------------------------
