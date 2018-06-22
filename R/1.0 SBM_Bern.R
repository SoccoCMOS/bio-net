# -----------------------------------------------------------------------------
# PROJECT:
#    BISE : Linking soil foodweb and soil functions dynamics under
#    agricultural practices (des)intensification and different input levels
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# STEP:
#    0.2    Building a Bernoulli-SBM on binary interactions 
#           and extracting adjacency matrix between groups
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# NOTES:
#  
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# LIBRARIES:
library(blockmodels)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Working directory
setwd("C:/1. Mike/5. RH/1. Stagiaire/18_Labouyrie Maeva/v2/bio-net/Data")
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# data load
# -----------------------------------------------------------------------------
adjacence_esp_bin <- read.csv2("adjacence_esp_bin.csv", h = T, sep =";")
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Building trophic groups through Stochastic Block Modelling
# -----------------------------------------------------------------------------    
my_model <- BM_bernoulli("SBM", adjacence_esp_bin,
                        verbosity=6,
                        plotting=character(0),
                        exploration_factor=1.2,
                        explore_max=100,
                        autosave = 'SBM_bern',
                        ncores=detectCores())
my_model$estimate()
which.max(my_model$ICL)
output.combined = readRDS("SBM_bern") 
output.combined$estimate()
Q = which.max(output.combined$ICL) # Optimal number of groups

# Best model
best.model = output.combined$memberships[[Q]] #Q
output.combined$plot_parameters(Q)

# Memberships         
SBM_member = as.data.frame(cbind(row.names(spec_adj), best.model$Z))
SBM_member$class <- NULL
for(i in 1:nrow(SBM_member)){
  SBM_member$class[i] <-  which.max(as.matrix(SBM_member[i, 2:(Q+1)])) #la premiere colonne correspond aux noms d'espece (Id)
}
SBM_member <- SBM_member[, c(1, ncol(SBM_member))]
colnames(SBM_member) = c("retained_tax", "class")  ##44 grp
# -----------------------------------------------------------------------------  


# -----------------------------------------------------------------------------
# Building SBM groups - adjacency matrix
# -----------------------------------------------------------------------------
M = output.combined$model_parameters
SBM_adj=M[[Q]]$pi #Q: qui est le max de groupes
colnames(SBM_adj) <- c(1:Q)
rownames(SBM_adj) <- colnames(SBM_adj)

#Ajout des familles et du code FW
SBM_member$fam <- NULL
SBM_member$codeFW <- NULL
for (i in 1:nrow(SBM_member)) {
  sp <- SBM_member$retained_tax[i]
  sp <- as.character(sp[1])
  fam <- unique(metabar$family_name[which(metabar$retained_tax == sp)])
  SBM_member$fam[i] <- as.character(fam[1])
  FW <- unique(metabar$codeFW[which(metabar$retained_tax == sp)])
  SBM_member$codeFW[i] <- as.character(FW[1])
}
write.csv2(SBM_member, "SBM_Bern_class.csv")
write.csv2(SBM_adj, "SBM_Bern_adj.csv")
# -----------------------------------------------------------------------------
