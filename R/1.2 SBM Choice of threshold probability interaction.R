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
library(igraph)
library(reshape)
library(ggplot2)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Working directory
setwd("../data")
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# data load
# -----------------------------------------------------------------------------

SBM_adj0 <- read.csv2("adjgroups.csv", h = T, sep =";")
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------

#créer réseau à partir de SBM edge list & de differents seuil sur p d'interactions dans adj matrix 
# sur le métaréseau


# keeping interactions with p > 0.01
p <- seq(0.005,0.995,by=0.005)

densities<-vector(mode="integer",length=length(p))
edges<-vector(mode="integer",length=length(p))
nodes<-vector(mode="integer",length=length(p))

SBM_adj <- SBM_adj0[,-1]

cpt=1
for (i in p){
  
  SBM_adj_filt <- as.data.frame(ifelse(SBM_adj > i , 1, 0))
  
  condition <- as.data.frame(rowSums(SBM_adj_filt))
  condition$colSums <- colSums(SBM_adj_filt)
  condition$cond <- rowSums(condition)
  
  SBM_adj_filt <- SBM_adj_filt[condition$cond  > 0, condition$cond > 0]
  
  
  # Creation du graphe
  g<-graph.adjacency(as.matrix(SBM_adj_filt))

  
  densities[cpt]=graph.density(g)
  edges[cpt]=len(graph.edgelist())
  edges[cpt]=E(g)$number
  
  cpt<-cpt+1
}

densities

