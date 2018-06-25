library(vegan)

SBM_clust <- read.csv("SBM_Poiss_class.csv", sep=";", header = T)
hierarchical_clust <- read.csv("groupsW-10.csv", sep=";", header=T)
spectral_clust <- read.csv("10_communities.csv", sep=";", header=T)

control <- as.matrix(table(SBM_clust$codeFW, SBM_clust$class))
plot(table(specnumber(control, MARGIN = 2)), ylab = "", las=1)
