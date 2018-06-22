# -----------------------------------------------------------------------------
# PROJECT:
#    BISE : Linking soil foodweb and soil functions dynamics under
#    agricultural practices (des)intensification and different input levels
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# STEP:
#    0.1   Comparing classification methods
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# NOTES:
#  
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# LIBRARIES:
library(blockmodels)
library(igraph)
library(Rmisc)
library(vegan)
library(betalink)
library(MASS)
library(reshape2)
library(ggplot2)
library(dplyr)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Working directory
setwd("J:/PROJET RECONSTRUCTION DE RESEAUX/TEST MATRICE TUTEUR/Script 21062018")
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Building trophic groups through Stochastic Block Modelling
# -----------------------------------------------------------------------------    
    my_model <- BM_gaussian("SBM", spec_adj,
                            verbosity=6,
                            plotting=character(0),
                            exploration_factor=1.2,
                            explore_max=100,
                            autosave = 'SBM_gauss',
                            ncores=detectCores())
    my_model$estimate()
    which.max(my_model$ICL)
    output.combined = readRDS("SBM_gauss") 
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
    SBM_adj=M[[Q]]$mu #Q: qui est le max de groupes
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
    
    # keeping interactions with p > 0.01
    SBM_edges <- melt(SBM_adj)
    SBM_edges <- SBM_edges[SBM_edges$value>0.01,]
    SBM_adj   <- dcast(SBM_edges, Var1~Var2, sum)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# SBM-groups metanetworks
# -----------------------------------------------------------------------------
#repr?sentation du m?tar?seau
eg <- read.csv("edge.csv", h = T, sep = ";")
eg2 <- merge(member1, eg[,1:2], "codeFW")
eg3 <- round(aggregate(eg2$TL, list(class = eg2$class),  mean), 2)
colnames(eg3)[2] <- "TL"

metag <- graph_from_data_frame(el5[,1:2], vertices = eg3, directed = TRUE)
lay<-matrix(nrow = Q, ncol=2)
lay[,1] <- runif(n = Q, min = 0.2, max = 0.8)
lay[,2]<- V(metag)$TL
V(metag)$name <- ""
plot.igraph(metag, layout = lay, vertex.label.font=2, vertex.size= 6,
            edge.width = 1, edge.lty = 1, edge.arrow.size=0)    
    
    
    
    
    
# -----------------------------------------------------------------------------
# Building SBM-groups networks
# -----------------------------------------------------------------------------    
    mb_short <- metabar[metabar$retained_tax %in% SBM_member$retained_tax, c(8, 18:273)]
    mb_short <- merge(SBM_member, mb_short, "retained_tax")
    mb_class <- aggregate(mb_short[, -c(1:3)], list(class = mb_short$class), sum)
    mb_class_t0 <- as.data.frame(t(mb_class))
    colnames(mb_class_t0) <- c(1:max(SBM_member$class))
    mb_class_t0 <- mb_class_t0[-1,]
    mb_class_t1 <- aggregate(mb_class_t0, by = list(point = paste(fac$Couples, fac$Echantillon, sep = "_")), FUN = sum)
    
    #A l'?chelle de la parcelle
    mb_list_t <- split(mb_class_t1[,-1], with(mb_class_t1[,-1], 
                                              rep(c("A+", "B+", "B-", "A-", "C+", "C-", "E+", "E-", "D+", "D-", 
                                                    "G-", "G+", "F-", "F+", "H+", "H-"), each = 4)))
    mb_list <- lapply(mb_list_t, FUN = t)
    
    #Cr?ation des r?seaux ? l'?chelle du plot
    mynetwork <- function(metabar){
      step1 <- ifelse(metabar>1,1,0)
      step2 <- step1[rowSums(step1)>0,]
      interactions <- SBM_adj[SBM_adj$Var1 %in% rownames(step2) & SBM_adj$Var2 %in% rownames(step2),]
      g <- graph_from_data_frame(interactions, directed = TRUE)
    }
    ntwk_list <- lapply(mb_list, FUN = mynetwork)
    
    #betadiversit? des r?seaux
    mybeta <- network_betadiversity(ntwk_list,  bf = B10)
    hist(mybeta$S)
    mybeta$S1 <- cut(mybeta$S, breaks = c(0.75,0.82, 0.9, 1),right = FALSE)
    ggplot(data = mybeta, aes(x=i, y=j)) + 
      geom_tile(aes( fill = ST), colour = "white")#+
    #scale_fill_brewer(palette = "PRGn")
    
    network_betaplot(ntwk_list[[4]], ntwk_list[[15]])
# -----------------------------------------------------------------------------
    
    
    
    #Cr?ation des r?seaux ? l'?chelle du point
    mb_list_t <- split(mb_class_t0, with(mb_class_t0, point))
    mb_list <- lapply(mb_list_t, FUN = t)
    mynetwork <- function(metabar){
      step1 <- ifelse(metabar>1,1,0)
      step2 <- rownames(step1)[step1[,1]>0]
      interactions <- el5[el5$Var1 %in% step2 & el5$Var2 %in% step2,]
      g <- graph_from_data_frame(interactions, directed = TRUE)
    }
    ntwk_list <- lapply(mb_list, FUN = mynetwork)
    
    #Calcul des indices de r?seau pour chaque point
    mynetwork_indices <- function(g){
      density <- edge_density(g, loops = FALSE)
      ntwkdegree <- degree(g, mode = "total")
      meandegree<- mean(ntwkdegree)
      nb_edge <- nrow(interactions)
      nb_node <- length(step2)
      obj <- list(density = density, meandegree = meandegree, nb_edge = nb_edge, nb_node = nb_node)
      return(obj)
    }
    ntwk_indic <- lapply(ntwk_list, FUN = mynetwork_indices)
    indices <- data.frame(matrix(unlist(ntwk_indic, use.names = TRUE), nrow=64, byrow=T))
    colnames(indices) <- c("density", "mean_degree", "nb_edge", "nb_node")
    indices$LECA_id <- point
    indices$Plot <- rep(plot, each = 4)
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
    
    