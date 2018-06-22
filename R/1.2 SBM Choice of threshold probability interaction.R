# -----------------------------------------------------------------------------
# PROJECT:
#    BISE : Linking soil foodweb and soil functions dynamics under
#    agricultural practices (des)intensification and different input levels
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# STEP:
#    1.2    Choice of probabiliy threshold in SBM adjacency matrices
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# NOTES:
#  
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# LIBRARIES:
library()
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Working directory
setwd("C:/1. Mike/5. RH/1. Stagiaire/18_Labouyrie Maeva/v2/bio-net/Data")
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# data load
# -----------------------------------------------------------------------------
SBM_Gauss_adj <- read.csv2("SBM_Gauss_adj.csv", h = T, sep =";")
SBM_Bern_adj <- read.csv2("SBM_Bern_adj.csv", h = T, sep =";")
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------

#créer réseau à partir de SBM edge list & de differents seuil sur p d'interactions dans adj matrix 
# sur le métaréseau


# keeping interactions with p > 0.01
SBM_edges <- melt(SBM_adj)
SBM_edges <- SBM_edges[SBM_edges$value>0.01,]
SBM_adj   <- dcast(SBM_edges, Var1~Var2, sum)

#représentation + test stat?

