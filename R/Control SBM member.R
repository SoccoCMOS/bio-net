library(vegan)

##### Visual control of memberships - expert analysis #######

groups <- read.csv("../data/input/Q_30_groups.csv", sep=";", header = T)
#hierarchical_clust <- read.csv("groupsW-10.csv", sep=";", header=T)
#spectral_clust <- read.csv("10_communities.csv", sep=";", header=T)
names <- read.csv2("../data/input/names.csv", header=T)

SBMcol1 <- names$x
SBMcol2 <- groups$V1
SBM <- data.frame(SBMcol1, SBMcol2)
colnames(SBM) <- c("retained_tax", "class")
SBM_clust <- SBM

#Ajout des familles et du code FW
metabar     <- read.csv("../data/input/metabar.csv", sep = ";", h = T) 

SBM_clust$fam <- NULL
SBM_clust$codeFW <- NULL
for (i in 1:nrow(SBM_clust)) {
  sp <- SBM_clust$retained_tax[i]
  sp <- as.character(sp[1])
  fam <- unique(metabar$family_name[which(metabar$retained_tax == sp)])
  SBM_clust$fam[i] <- as.character(fam[1])
  FW <- unique(metabar$codeFW[which(metabar$retained_tax == sp)])
  SBM_clust$codeFW[i] <- as.character(FW[1])
}

control <- as.matrix(table(SBM_clust$codeFW, SBM_clust$class))
plot(table(specnumber(control, MARGIN = 2)), ylab = "", las=1)
table(specnumber(control, MARGIN = 2))


###### Detail: qu'est-ce que fonctionne moins bien ? #####
ctr<-SBM_clust[,c("codeFW","class")]
detail=aggregate(ctr$codeFW,by=list(ctr$class),function(x) length(unique(x)))
