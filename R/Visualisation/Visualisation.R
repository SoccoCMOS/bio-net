###################################
#### VISUALISATION DES DONNEES ####
###################################

setwd(".") 

library(igraph)
library(reshape2)
library(vegan)
library(betalink)
library(ggplot2)

# -----------------------------------------------------------------------------
# Networks
# -----------------------------------------------------------------------------

###
# METARESEAU
###

###NIVEAU FAMILLE


###NIVEAU GENRE


###NIVEAU ESPECE
adjacence_esp_bin <- read.csv2("adjacence_esp_bin.csv", header=T, row.names=1)
colnames(adjacence_esp_bin) <- rownames(adjacence_esp_bin) ##Uniform row and col names
setdiff(colnames(adjacence_esp_bin), rownames(adjacence_esp_bin))
adjacence_esp_bin <- as.matrix(adjacence_esp_bin) ## Format accepted by igraph

metareseau <- graph.adjacency(adjacence_esp_bin, mode = "directed", weighted = NULL, diag = FALSE)
#metareseau <- simplify(metareseau, remove.multiple = F, remove.loops = T)
plot(metareseau, edge.arrow.size=.2, edge.curved =.1,
     vertex.size = 3, vertex.label=NA)

#Metrics 
connectivity <- edge_density(metareseau, loops = FALSE) #connectivity = number of observed links divided by the number of possible links without cannibalism (L/S(S-1))
ntwkdegree <- degree(metareseau, mode = "total") #number of adjacent edge of a vertex
ntwkindegree <- degree(metareseau, mode = "in") #number of adjacent edge of a vertex
ntwkoutdegree <- degree(metareseau, mode = "out") #number of adjacent edge of a vertex
meandegree<- mean(ntwkdegree) 
meanindegree<- mean(ntwkindegree) #linkage density
meanoutdegree<- mean(ntwkoutdegree)
nb_edge <- gsize(metareseau) #number of edges
nb_node <- gorder(metareseau) #number of nodes or vertices



filt=c("symbiotic")

###
# RESEAU EXPERT
###
eel <- read.csv2("expert_edge_list.csv")
el <- subset(eel,!(eel$type %in% filt))
el <- as.matrix(el[,c("TG1", "TG2")])

expg<-graph_from_edgelist(el,directed=TRUE)

eg_connectivity <- edge_density(expg, loops = FALSE) #linkage density
eg_ntwkdegree <- degree(expg, mode = "total") #number of adjacent edge of a vertex
eg_ntwkindegree <- degree(expg, mode = "in") #number of adjacent edge of a vertex
eg_ntwkoutdegree <- degree(expg, mode = "out") #number of adjacent edge of a vertex
eg_meandegree<- mean(eg_ntwkdegree)
eg_meanindegree<- mean(eg_ntwkindegree) #linkage density
eg_meanoutdegree<- mean(eg_ntwkoutdegree)
eg_nb_edge <- gsize(expg) #number of edges
eg_nb_node <- gorder(expg) #number of nodes or vertices



###
# RESEAU SBM
###
# A partir de la matrice d'adjacence des "31" groupes

#Lecture fichier SBM comportant la matrice d'adjacence entre groupes à un threshold de ??
matadjSBM <- read.csv("../filtered_sbm_adj_mat.csv", dec=",", sep=";", header=T, row.names = 1)

colnames(matadjSBM) <- rownames(matadjSBM)
matadjSBM <- as.matrix(matadjSBM)


# keeping interactions with p > 0.01 ### SCRIPT P-THRESHOLD
# p=0.01
# SBM_edges_in <- melt(matadjSBM)
# SBM_edges_in$value<-ifelse(SBM_edges_in$value>p,1,0)
# 
# SBM_edges <- subset(SBM_edges_in, SBM_edges_in$value==1)
# 
# SBM_adj_bin   <- as.matrix(dcast(SBM_edges_in, Var1~Var2, fill=0)[,-1]) 

SBM_adj_bin <- matadjSBM
reseau_sbm <- graph.adjacency(SBM_adj_bin, mode = "directed", weighted = NULL, diag = FALSE)

sbm_connectivity <- edge_density(reseau_sbm, loops = FALSE) #linkage density
sbm_ntwkdegree <- degree(reseau_sbm, mode = "total") #number of adjacent edge of a vertex
sbm_ntwkindegree <- degree(reseau_sbm, mode = "in") #number of adjacent edge of a vertex
sbm_ntwkoutdegree <- degree(reseau_sbm, mode = "out") #number of adjacent edge of a vertex
sbm_meandegree<- mean(sbm_ntwkdegree) 
sbm_meanindegree<- mean(sbm_ntwkindegree) #linkage density
sbm_meanoutdegree<- mean(sbm_ntwkoutdegree)
sbm_nb_edge <- gsize(reseau_sbm) #number of edges
sbm_nb_node <- gorder(reseau_sbm) #number of nodes or vertices

#reseau_SBM <- simplify(reseau_SBM, remove.multiple = F, remove.loops = T)
plot(reseau_sbm, edge.arrow.size=.2, edge.curved =.1,
     vertex.size = 3, vertex.label=NA)


SBM_member <- read.csv2("C:/Users/simoussi/ownCloud/Ma_these/bio-net/data/input/label_SBM30.csv", header=T)
SBM_member <- SBM_member[,c("retained_tax","class", "fam", "codeFW")]
metabar <- read.csv2("C:/Users/simoussi/ownCloud/Ma_these/bio-net/data/input/metabar.csv", header=T)
orig_fac <- read.csv2("C:/Users/simoussi/ownCloud/Ma_these/bio-net/data/input/fac.csv", header = T)

# -----------------------------------------------------------------------------
# Building SBM-groups networks
# -----------------------------------------------------------------------------    
mb_short <- metabar[metabar$retained_tax %in% SBM_member$retained_tax, c(8, 18:273)]
mb_long <- merge(SBM_member, mb_short, "retained_tax")
mb_class <- aggregate(mb_long[, -c(1:4)], list(class = mb_long$class), sum)
mb_class_tf <- as.data.frame(t(mb_class))
colnames(mb_class_tf) <- c(0:max(SBM_member$class))
mb_class_t0 <- mb_class_tf[-1,]

### mb_class_t1 [i,j] = Number of reads of MOTU classified in class j seen in point i => group per point (replicates up to point)
fac<-orig_fac[order(orig_fac$Couples),]
mb_class_t1 <- aggregate(mb_class_t0, by = list(point = paste(fac$Couples, fac$Echantillon, sep = "_")), FUN = sum)

#A l'?chelle de la parcelle --> grouper les 4 échantillons à la parcelle (points up to plots)
mb_list_t <- split(mb_class_t1[,-1], with(mb_class_t1[,-1], 
                                          rep(c("A-", "A+", "B-", "B+", "C-", "C+", "D+", "D-", "E-", "E+", 
                                                "F-", "F+", "G-", "G+", "H-", "H+"), each = 4)))
mb_list <- lapply(mb_list_t, FUN = t)

#Cr?ation des r?seaux ? l'?chelle du plot
mynetwork <- function(mat){
  step1 <- ifelse(mat>1,1,0) ## Binarize, if at least one read of a taxa in group on row, say true else false ???? TODO: >1 or 100 or 0 ?
  step2 <- step1[rowSums(step1)>0,] ## Number of points where the considered class is read
  interactions <- SBM_edges[SBM_edges$Var1 %in% rownames(step2) & SBM_edges$Var2 %in% rownames(step2),]
  g <- graph_from_data_frame(interactions, directed = TRUE)
  g <- simplify(g, remove.multiple = F, remove.loops = T)
}


ntwk_list_plot <- lapply(mb_list, FUN = mynetwork)

#betadiversit? des r?seaux
mybeta <- network_betadiversity(ntwk_list_plot,  bf = B11)
hist(mybeta$S)
mybeta$S1 <- cut(mybeta$S, breaks = c(0.75,0.82, 0.9, 1),right = FALSE)
ggplot(data = mybeta, aes(x=i, y=j)) + 
  geom_tile(aes( fill = ST), colour = "white")#+
#scale_fill_brewer(palette = "PRGn")

network_betaplot(ntwk_list_plot[[1]], ntwk_list_plot[[12]]) ### Vert première parcelle (réseau) exclusivement, bleu 2ème parcelle, gris les deux
# -----------------------------------------------------------------------------



#Cr?ation des r?seaux ? l'?chelle du point
mb_list_t <- split(mb_class_t1, with(mb_class_t1, point))
mb_list_point <- lapply(mb_list_t, FUN = t)
ntwk_list_point <- lapply(mb_list_point, FUN = mynetwork)

#Calcul des indices de r?seau pour chaque point+plot
mynetwork_indices <- function(g){
  connectivity <- edge_density(g, loops = FALSE)
  ntwkdegree <- degree(g, mode = "total")
  meandegree<- mean(ntwkdegree)
  nb_edge <- ecount(g)
  nb_node <- vcount(g)
  density <- nb_edge/nb_node
  obj <- list(density= density, connectivity = connectivity, meandegree = meandegree, nb_edge = nb_edge, nb_node = nb_node)
  return(obj)
}

ntwk_indic_point <- lapply(ntwk_list_point, FUN = mynetwork_indices)
ntwk_indic_plot <- lapply(ntwk_list, FUN = mynetwork_indices)

indices_point <- data.frame(matrix(unlist(ntwk_indic_point, use.names = TRUE), nrow=64, byrow=T))
colnames(indices_point) <- c("density","connectivity", "mean_degree", "nb_edge", "nb_node")
indices_point$LECA_id <- point
indices_point$Plot <- rep(plot, each = 4)
fac <- read.csv("fac.csv", h = T, sep = ";")
ind <- merge(indices, fac,  "LECA_id")
density <- summarySE(data = ind, measurevar = "density", groupvars = "Couples", na.rm = T)
nb_edge <- summarySE(data = ind, measurevar = "nb_edge", groupvars = "Couples", na.rm = T)
nb_node <- summarySE(data = ind, measurevar = "nb_node", groupvars = "Couples", na.rm = T)
mean_degree <- summarySE(data = ind, measurevar = "mean_degree", groupvars = "Plot", na.rm = T)


#Repr?sentation des indices de r?seau
ggplot(density, aes(x = Couples, y = density))+
  geom_point()+
  geom_errorbar(aes(ymin = density - se, ymax = density + se))+
  geom_vline(xintercept = c(2.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5), size = 1.2, 
             linetype = 2)+
  geom_vline(xintercept = c(4.5, 10.5), size = 1.5, linetype = 1, colour = "red")+
  ylab("Density")+
  xlab("")+
  theme(panel.background = element_rect(fill = 'white', colour = 'red'))

ggplot(mean_degree, aes(x = Plot, y = mean_degree))+
  geom_point()+
  geom_errorbar(aes(ymin = mean_degree - se, ymax = mean_degree + se))+
  geom_vline(xintercept = c(2.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5), size = 1.2, 
             linetype = 2)+
  geom_vline(xintercept = c(4.5, 10.5), size = 1.5, linetype = 1, colour = "red")+
  ylab("mean_degree")+
  xlab("")+
  theme(panel.background = element_rect(fill = 'white', colour = 'red'))




#######
# Fonctions microbiennes mesur?es
#######
funct <- read.csv(file = "functCN.csv", sep = ";")
fun <- aggregate(funct[,-c(1:4)], list(funct$couple), mean)
ind_agg <- aggregate(ind[,2:3], list(ind$Couples), mean)


fun1 <- summarySE(data = funct, measurevar = "C.sol", groupvars = "couple", na.rm = T)

site <- c(rep("MONS", 4), rep("LUS", 6), rep("THX", 6))
plot(ratio_dec_min$ratio_dec_min~fun1$C.sol, 
     col = c(rep("darkblue", 4), rep("darkred", 6), rep("darkgreen", 6)),  
     pch = 3, xlab = "Soil C content (mgC kg soil)", ylab = "decomposer:miner ratio")


plot(fun1$C.sol~ratio_dec_min$ratio_dec_min, 
     col = c(rep("darkblue", 4), rep("darkred", 6), rep("darkgreen", 6)), 
     pch = 3, ylab = "C mineralization speed (mgC d-1 kg soil)", xlab = "decomposer:miner ratio")

BEF <- cbind (ind_agg, fun)
cormat <- round(cor(BEF[,-c(1, 4)]),2)
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="Pearson\nCorrelation") +
  #theme_minimal()+ 
  #theme(axis.text.x = element_text(angle = 45, vjust = 1, 
  #                                 size = 12, hjust = 1))+
  coord_fixed()

