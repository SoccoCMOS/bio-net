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

SBM_adj <- read.csv2("../data/Q_33_groups.csv", h = T, sep =";", row.names = 1)
rownames(SBM_adj)=colnames(SBM_adj)
# -----------------------------------------------------------------------------
##### Règle de décision sur le choix de p ######
## Criterion 1: Keep alpha ratio of nodes (trophic groups)
th_n=1 ## Max tolerated loss ratio of nodes
max_nbintra=0 ## Max accepted number of intraguild interactions
## Criterion 2: Keep graph connexity
## Criterion 3: No intraguild interactions
# -----------------------------------------------------------------------------

#créer réseau à partir de SBM edge list & de differents seuil sur p d'interactions dans adj matrix 
# sur le métaréseau

# keeping interactions with p > 0.01
p_range <- seq(0,0.995,by=0.001)

# Following metrics are computed to keep track of changes in graph entropy as we filter edges based on weight threshold p
densities<-vector(mode="integer",length=length(p_range)) # Graph density = nb_edges/nb_nodes
# edges<-vector(mode="integer",length=length(p)) # Number of edges (hear interactions)
# nodes<-vector(mode="integer",length=length(p)) # Number of connected nodes

n=length(SBM_adj)

optim_reached=FALSE
node_preserv=TRUE #Criterion 1
intrag=TRUE # criterion 3
connex=TRUE # criterion 2
cpt=1
while (!optim_reached){
  cpt=cpt+1
  p=p_range[cpt]
  ## Binarize
  SBM_adj_filt <- as.data.frame(ifelse(SBM_adj > p , 1, 0))
  ## Filter
  condition <- as.data.frame(rowSums(SBM_adj_filt))
  condition$colSums <- colSums(SBM_adj_filt)
  condition$cond <- rowSums(condition)
  SBM_adj_filt <- SBM_adj_filt[condition$cond  > 0, condition$cond > 0]
  
  # Creation du graphe
  g<-graph.adjacency(as.matrix(SBM_adj_filt))

  # Re-compute metrics
  densities[cpt-1]=graph.density(g)
  #nodes[cpt]=vcount(g)
  #edges[cpt]=ecount(g)
  
  # Node preservation criterion (1)
  node_preserv=(vcount(g)/n)>=th_n 

  # Check for intraguild criterion (3)
  intra=lapply(1:length(SBM_adj_filt),FUN=function(x) SBM_adj[x,x]>p)
  n_intrag=length(which(intra==TRUE))
  intrag=n_intrag>max_nbintra
  
  # Check for connexity criterion (2)
  clu<-components(g)
  connex=clu$no==1
  
  if(!(node_preserv & intrag & connex)){
    optim_reached=TRUE
    p=p_range[cpt-1]
  }
}

densities

par(mfrow=c(2,2))
plot(densities~p_range)

# plot(nodes~p_range)
# plot(edges~p_range)
# edges2 <- edges[2:length(edges)]-edges[1:(length(edges)-1)]
# plot(edges2~p_range[-197])

#########################################################

SBM_adj_filt <- as.data.frame(ifelse(SBM_adj > p , 1, 0))

condition <- as.data.frame(rowSums(SBM_adj_filt))
condition$colSums <- colSums(SBM_adj_filt)
condition$cond <- rowSums(condition)

SBM_adj_filt <- SBM_adj_filt[condition$cond  > 0, condition$cond > 0]
write.csv2(SBM_adj_filt,"filtered_sbm_adj_mat.csv")
