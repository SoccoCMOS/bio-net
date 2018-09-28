# -----------------------------------------------------------------------------
# PROJECT:
#    BISE : Linking soil foodweb and soil functions dynamics under
#    agricultural practices (des)intensification and different input levels
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# STEP:
#    5.0   Analyzing metanetworks
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# NOTES:
# 
install.packages("igraph")
install.packages("RNewsflow")
install.packages("Rmisc")
install.packages("betalink")

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# LIBRARIES:
library(igraph)
library(RNewsflow)
library(reshape2)
library(Rmisc)
library(betalink)
library(ggplot2)
library(dplyr)
library(stringr)

# -----------------------------------------------------------------------------
################################################################# FUNCTIONS  ###################################################################
#################################################### BEGIN(Brown/green pathways) #############################################################

get_brown=function(clust,meth){
  if(meth=="expert"){
    roots=resources
  }
  if(meth=="none"){
    roots=resources
  }
  if(meth=="sbm"){
    roots_unform=unique(subset(clust,clust$retained_tax %in% resources)$class)
    roots=unlist(lapply(roots_unform,function(x) paste("V",x,sep="")))
  }
  
  return(roots)
}

get_green=function(clust,meth){
  green_roots_df=clust[grep("Phototroph::*",x=clust$codeFW),]
  
  if(meth=="expert"){
    roots=green_roots_df$codeFW
  }
  if(meth=="none"){
    roots=green_roots_df$retained_tax
  }
  if(meth=="sbm"){
    roots=unlist(lapply(unique(green_roots_df$class),function(x) paste("V",x,sep="")))
  }
  
  return(roots)
}


pathways_stat<-function(net,sources){    ### Computes pathways statistics => provide phototroph on sources for green pathways (resources for brown respectively)
  if(vcount(net)==0){
    print("Warning empty or disconnected network")
    brg_stat=list(sum_pathsize=0,mean_pathsize=0,median_pathsize=0,sd_pathsize=0,covered_taxa=NULL,nbcovered_taxa=0)
  }
  else{
    node_sources<-V(net)[name %in% sources]
    deb <<- net
    sums=vector(mode="numeric",length=length(node_sources))
    moys=vector(mode="numeric",length=length(node_sources))
    meds=vector(mode="numeric",length=length(node_sources))
    sds=vector(mode="numeric",length=length(node_sources))
    maxs=vector(mode="numeric",length=length(node_sources))
    
    ncov=vector(mode="character")
    
    cpt=0
    for (root in node_sources){
      cpt=cpt+1
      d=bfs(net, root, neimode = c("in"), 
            unreachable = FALSE, restricted = NULL, order = TRUE, rank = FALSE,
            father = FALSE, pred = FALSE, succ = FALSE, dist = TRUE,
            callback = NULL, extra = NULL, rho = parent.frame()) 
      
      roads=d$dist[which(!is.na(d$dist))]
      sums[cpt]=sum(roads)
      moys[cpt]=mean(roads)
      meds[cpt]=median(roads)
      sds[cpt]=sd(roads)
      maxs[cpt]=max(roads)
      ncov=union(ncov,names(which(!is.na(d$dist))))
    }
    
    nbcov=length(ncov)
    brg_stat=list(sum_pathsize=sum(sums),mean_pathsize=mean(moys),median_pathsize=median(meds),sd_pathsize=sd(sds),covered_taxa=ncov,nbcovered_taxa=nbcov,max_dist=max(maxs))
  }
  return(brg_stat)
}

ratio<-function(a,b){
  if(is.na(a) || is.na(b)){
    rat=0
  }
  else{
    if(b==0){
      rat=0
    }
    else{
      rat=a/b
    }
  }
}

brown_green <- function(net,green,brown){   ### Analysis of a single network
  B=pathways_stat(net,brown)
  G=pathways_stat(net,green)
  
  ### Brown/Green statistics
  comb=data.frame(source=c("brown/green"),sum=ratio(B$sum_pathsize,G$sum_pathsize), med=ratio(B$median_pathsize,G$median_pathsize),mean=ratio(B$mean_pathsize,G$mean_pathsize),sd=ratio(B$sd_pathsize,G$sd_pathsize),nbcovered=ratio(B$nbcovered_taxa,G$nbcovered_taxa))
}

brown_green_detail <- function(net,green,brown){   ### Analysis of a single network
  B=pathways_stat(net,brown)
  G=pathways_stat(net,green)
  
  ### Brown/Green statistics
  BR=data.frame(source=c("brown"),sum=B$sum_pathsize, med=B$median_pathsize,mean=B$mean_pathsize,sd=B$sd_pathsize,nbcovered=B$nbcovered_taxa)
  GR=data.frame(source=c("green"),sum=G$sum_pathsize, med=G$median_pathsize,mean=G$mean_pathsize,sd=G$sd_pathsize,nbcovered=G$nbcovered_taxa)
  out=list(BR=BR,GR=GR)
}
#################################################### END (Brown/green pathways computing) #############################################################
######################################################### Params ########################################################

############### INPUT FILES and WORKING DIRECTORIES #################
setwd("../data/input/taxo_metabar/adj")
net_folders=c("sbm/t","expert","none")
clust_folder=c("../tax_sbm_exp")
files=c("family.csv","genus.csv","species.csv")
mb_folder="../mb"
fac_file="../../fac.csv"

out_folder=c("../out_alpha/sbm","../out_alpha/expert","../out_alpha/none")
files_point=c("point_family.csv","point_genus.csv","point_species.csv")

## RÃ©partition des parcelles, plots, Ã©chantillons, rÃ©plicats
orig_fac=read.csv2(fac_file)
fac<-orig_fac#[order(orig_fac$Couples),]

### CONSTANTS USED ###
colu=c("class","codeFW","retained_tax")
method=c("sbm","expert","none")
taxo=c("family","genus","species")
resources=c("Resource::Decaying_wood","Resource::Dung","Resource::MOS","Resource::Litter")

#### For beta div
indice=B10  #Distance metric Jaccard

### Initializations
brown_roots=list(length=9)
green_roots=list(length=9)
deb=NULL
cpt=0
results=list()


for (j in 1:(length(colu))){  ### Clustering method variation
  cl=colu[j]
  
  for (i in 1:length(taxo)){ ### Taxonomic scale variation
  # -----------------------------------------------------------------------------
  # data load
  # -----------------------------------------------------------------------------
  ## Matrice d'adjacence du rÃ©seau
   net=read.csv2(paste(net_folders[j],files[i],sep="/"),row.names = 1)
  ## Matrice metabarcoding projetÃ©e au niveau taxonomique i
   mb=read.csv2(paste(mb_folder,files[i],sep="/"),row.names = 1)
  ## Labels du clustering
   clust=read.csv2(paste(clust_folder,files[i],sep="/"),row.names = 1)
  # -----------------------------------------------------------------------------
  ### Edge list
  ## Uniformiser les noms des lignes et colonnes (formatage)
  #row.names(net)=as.character(1:length(net)) 
  colnames(net)=row.names(net)
  ## CrÃ©ation de la edge list Ã  partir de la matrice d'adjacence
  net_edges=melt(as.matrix(net))
   
  ### Projection du m?tar?seau sur la parcelle ###
  ## Filtrer les taxons qui interagissent uniquement + leurs relevÃ©s de read metabar 
  mb_short <- mb[mb$retained_tax %in% clust$retained_tax, c(273, 17:272)]
  ## Ajouter les informations de clustering: classe attribuÃ©e (codeFW pour expert,retained_tax pour none,class pour SBM)
  ## => Dans mb_long: pour chaque taxon (au niveau i), on aura les reads dans les 256 points d'Ã©chantillonnages (16*4*4) + diffÃ©rentes classifications
  mb_long <- merge(clust, mb_short, "retained_tax")
  ## Aggregate number of reads by summing the value of all rows with same class/codeFW/retained_taxa
  ## mb_class = aggregated metabar reads per class
  mb_class <- aggregate(mb_long[, -c(1:3)], list(class = mb_long[,cl]), sum)
  ## transpose to have metabar sampling points on rows and classes on cols
  mb_class_t0 <- as.data.frame(t(mb_class[,-1]))
  # ## Supprimer premiÃ¨re ligne (ancienne colonne index)
  # mb_class_t0<-mb_class_tf[-1,]
  
  ## Rename columns
  colnames(mb_class_t0) <- colnames(net)  
  
  ### mb_class_t1 [i,j] = Number of reads of MOTU classified in class j seen in point i => group per point (replicates up to point)
  ## Pour chaque plot (parcelle + traitement)
  mb_class_t1 <- aggregate(mb_class_t0, by = list(point = paste(fac$Couples, fac$Echantillon, sep = "_")), FUN = sum)
  
  #A l'echelle de la parcelle --> grouper les 4 ?chantillons ? la parcelle (points up to plots)
  mb_list_t <- split(mb_class_t1[,-1], with(mb_class_t1,point))
  mb_list <- lapply(mb_list_t, FUN = t)
  
  mb_list_t_plot <- split(mb_class_t1[,-1], with(mb_class_t1[,-1], 
                                            rep(c("A+", "B+", "B-", "A-", "C+", "C-", "E+", "E-", "D+", "D-", 
                                                  "G-", "G+", "F-", "F+", "H+", "H-"), each = 4)))
  mb_list_plot <- lapply(mb_list_t_plot, FUN = t)
  
  #Creation des reseaux a l'echelle du plot
  mynetwork <- function(mat){
    step1 <- ifelse(mat>1,1,0) ## Binarize, if at least one read of a taxa in group on row, say true else false ???? TODO: >1 or 100 or 0 ?
    step2 <- as.data.frame(step1[rowSums(step1)>0,]) ## Number of points where the considered class is read
    interactions <- net_edges[net_edges[,1] %in% rownames(step2) & net_edges[,2] %in% rownames(step2),]
    filt_inter <- subset(interactions,interactions$value>0)
    if(length(filt_inter)==0){
      print("stop")
    }
    g <- graph_from_data_frame(filt_inter, directed = TRUE)
    g <- simplify(g, remove.multiple = F, remove.loops = T)
  }
  
  ntwk_list_point <- lapply(mb_list, FUN = mynetwork)
  ntwk_list_plot <- lapply(mb_list_plot, FUN = mynetwork)
  
   
  ### Calcul des indices pour une analyse en Alpha ###
  mynet_ind <- function(g){
    # ntwkdegree <- degree(g, mode = "total")
    # meandegree<- mean(ntwkdegree)
    nb_edge <- ecount(g)
    nb_node <- vcount(g)
    connectivity <- nb_edge/(nb_node * (nb_node-1))
    density <- nb_edge/nb_node
    obj <- list(connectivity=connectivity, density= density, nb_edge = nb_edge, nb_node = nb_node)
    return(obj)
  }
  
  ntwk_indic_point <- lapply(ntwk_list_point, FUN = mynet_ind)
  indices_point <- data.frame(matrix(unlist(ntwk_indic_point, use.names = TRUE), nrow=64, byrow=T))
  colnames(indices_point) <- c("connectivity","density", "nb_edge", "nb_node")
  
  ntwk_indic_plot <- lapply(ntwk_list_plot, FUN = mynet_ind)
  indices_plot <- data.frame(matrix(unlist(ntwk_indic_plot, use.names = TRUE), nrow=16, byrow=T))
  colnames(indices_plot) <- c("connectivity","density", "nb_edge", "nb_node")  
  
  #### Current results
  res=list()
  
  res$nets_plot=ntwk_list_plot
  res$nets_point=ntwk_list_point
  
  #### Analyse en Alpha
  res$meth=net_folders[j]
  res$tax=taxo[i]
  res$metpoints=indices_point
  res$metplots=indices_plot
  
  ## Get brown and green roots
  brown_roots[[(j-1)*3+i]]=get_brown(clust,meth=method[j]) ###clust différencie family/genus/species
  green_roots[[(j-1)*3+i]]=get_green(clust,meth=method[j]) ###clust différencie family/genus/species
  
  ## Brown_green pathways stats
  res$brg_plot$brown=lapply(ntwk_list_plot,function(x) brown_green_detail(x,green_roots[[(j-1)*3+i]],brown_roots[[(j-1)*3+i]])$BR)
  res$brg_plot$green=lapply(ntwk_list_plot,function(x) brown_green_detail(x,green_roots[[(j-1)*3+i]],brown_roots[[(j-1)*3+i]])$GR)
  res$brg_point=lapply(ntwk_list_point,function(x) brown_green(x,green_roots[[(j-1)*3+i]],brown_roots[[(j-1)*3+i]]))
  
  #### Analyse en Beta
  #res$beta_plot=network_betadiversity(ntwk_list_plot,  bf = indice)
  #res$beta_point=network_betadiversity(ntwk_list_point,  bf = indice)
  
  results[[(j-1)*3+i]]=res
  
  write.csv2(indices_point,paste(out_folder[j],files_point[i],sep="/"))
  write.csv2(indices_plot,paste(out_folder[j],files[i],sep="/"))
  }
}



########################################### BEGIN (Comparative analysis computing and plotting functions) ##############################################################################
lin_mod=function(X,Y){    ##Linear model of Y~X   ### Fonction qui retourne la p-value, le R² de Y en fonction de X
  p=cor.test(X,Y)$p.value
  rsq=(cor.test(X,Y)$estimate)^2
  stat=list(R=rsq,pvalue=p)
}


### Plot Y~X ###
plot_relative <- function(X,Y,xl="",yl="",t=""){   ### Plot Y~X
  chart=plot(Y~X,xlab=xl,ylab=yl,main=t)
  print(chart)
}


###### Alpha plots #####
alpha_plot=function(x,y,taxonomic_level,titre){   ### BarPlot le R² par méthode de clustering et par niveau taxonomique
  # Grouped Bar Plot
  data <- data.frame(x,y,taxonomic_level)
  
  colnames(data)=c("Clustering method","R²","Taxonomic level")
  
  g=ggplot(data, aes(factor(x), y, fill = taxonomic_level)) + 
    geom_bar(stat="identity", position = "dodge") + 
    scale_fill_brewer(palette = "Set1") +
    labs(x="Clustering method", y="R²", title=titre)
  
  print(g)
}

#### Smooth graphs ####
alpha_smooth=function(x,y,z,tr,s){ ### Graphe de Y en fonction de X pour Z (méthode) et par site (shapes)
  data=data.frame(x,y,z)
  g=ggplot(data,aes(x=x,y=y,color=z,shape=s))+
    geom_point()+
    geom_abline(slope = 1, intercept=0)+
    geom_smooth(method=lm,se=TRUE,fullrange=FALSE) +
    labs(title=tr,x=paste("Level 1",tr,sep = " "),y=paste("Level 2",tr,sep=" "))
  
  print(g)
}

#lmetrics= list of 9 vectors containing the values of the metric for all combinations of method and taxo
global_analysis <- function(m="metric_name",s="scale",lmetrics,sh){ 
  
  #### R Square plot
  sgf=lin_mod(lmetrics[[1]],lmetrics[[2]])  ##R² entre family et genus pour SBM
  egf=lin_mod(lmetrics[[4]],lmetrics[[5]])  ## pour Expert
  ngf=lin_mod(lmetrics[[7]],lmetrics[[8]])  ## pour None
  
  ssg=lin_mod(lmetrics[[2]],lmetrics[[3]])  ##R² entre genus et species pour SBM
  esg=lin_mod(lmetrics[[5]],lmetrics[[6]])
  nsg=lin_mod(lmetrics[[8]],lmetrics[[9]])
  
  y=c(ssg$R,esg$R,nsg$R,sgf$R,egf$R,ngf$R) ## vecteurs des R² pour les 9 combi
  x=c(rep(c("SBM","EXPERT","NONE"),2))  ### Labels des méthodes
  z=c(rep("Species_Genus",3), rep("Genus_Family",3))  ### Labels des couples de taxo à comparer
  
  alpha_plot(x,y,z,paste(m,"at",s,"level",sep=" ")) ### Barplot
  
  #### Relative plot => plot l'une en fonction de l'autre
  plot_relative(lmetrics[[3]],lmetrics[[2]],paste(m,"at",s,"level","SBM","Species"),paste(m,"at",s,"level","SBM","Genus"))
  plot_relative(lmetrics[[2]],lmetrics[[1]],paste(m,"at",s,"level","SBM","Genus"),paste(m,"at",s,"level","SBM","Family"))
  
  plot_relative(lmetrics[[6]],lmetrics[[5]],paste(m,"at",s,"level","EXPERT","Species"),paste(m,"at",s,"level","EXPERT","Genus"))
  plot_relative(lmetrics[[5]],lmetrics[[4]],paste(m,"at",s,"level","EXPERT","Genus"),paste(m,"at",s,"level","EXPERT","Family"))
  
  plot_relative(lmetrics[[9]],lmetrics[[8]],paste(m,"at",s,"level","NONE","Species"),paste(m,"at",s,"level","NONE","Genus"))
  plot_relative(lmetrics[[8]],lmetrics[[7]],paste(m,"at",s,"level","NONE","Genus"),paste(m,"at",s,"level","NONE","Family"))
  
  #### Smoothed graph - species/genus
  sx=c(lmetrics[[3]],lmetrics[[6]],lmetrics[[9]]) ###Valeurs de la métrique pour species
  sz=rep(c("SBM","EXPERT","NONE"),each=length(lmetrics[[3]])) ### Méthode
  sy=c(lmetrics[[2]],lmetrics[[5]],lmetrics[[8]]) ###Valeurs de la métrique pour genus
  alpha_smooth(sx,sy,sz,m,sh)
  
  #### Smoothed graph - genus/family
  sx2=c(lmetrics[[2]],lmetrics[[5]],lmetrics[[8]])  ###Valeurs de la métrique pour genus
  sz2=rep(c("SBM","EXPERT","NONE"),each=length(lmetrics[[2]])) ### Méthode
  sy2=c(lmetrics[[1]],lmetrics[[4]],lmetrics[[7]])  ###Valeurs de la métrique pour family
  alpha_smooth(sx2,sy2,sz2,m,sh)
}

brg_ratios<-function(ld){ #l is a list of brown green metrics for a set of plots ==> mettre les brg pour tous les plots/points dans un seul dataframe
  mrg=ldply(ld, data.frame)
  out=mrg[,-c(1,2)] 
}
########################################### END (ALPHA ANALYSIS functions) ##############################################################################
##I. Création/Initialization des variables
######################################################### Alpha ##############################################################################
sh_pn=rep(c(rep("Mons",16),rep("Lusignan",24),rep("Theix",24)),3) ### Infos sur les parcelles (nom)
sh_pl=rep(c(rep("Mons",4),rep("Lusignan",6),rep("Theix",6)),3) ### //

#### (métrique= c pour connectivité, d pour densité, e pour edges, n pour nodes, sumbg pour ratio de la somme des chemins brown/green ....)(liste=l)(échelle=pn pour point / pl pour plot)
### Dans chaque liste, il y a 9 vecteurs, chaque vecteur correspond à une combinaison méthode+taxo
### Pour chaque métrique, deux listes, une par échelle (plot/point)

clpn=list(length=length(results)) 
clpl=list(length=length(results))

dlpn=list(length=length(results))
dlpl=list(length=length(results))

elpn=list(length=length(results))
elpl=list(length=length(results))

nlpn=list(length=length(results))
nlpl=list(length=length(results))

sumbglpn=list(length=length(results))  ##Ratio brown/green de la Somme des longueurs des chemins brown green en partant de toutes les sources
sumbglpl=list(length=length(results))

sdbglpn=list(length=length(results))  ##Ratio de l'Ecart type de la longueur des chemins
sdbglpl=list(length=length(results))

mnbglpn=list(length=length(results)) ##Ratio de la Moyenne de la longueur des chemins
mnbglpl=list(length=length(results))

covsbglpn=list(length=length(results))  ##Ratio du Nombre d'espèces atteignables 
covsbglpl=list(length=length(results))

maxsbglpn=list(length=length(results)) ##Ratio de la longueur maximale du chemin brown/green
maxsbglpl=list(length=length(results))

######################################################### Beta ##############################################################################

bpl=list(length=length(results))  ## béta à l'échelle du plot (liste de taille 9, 9 vecteurs de béta entre pairs des 16 plots)
#bpn=list(length=length(results)) ## idem échelle point (64 points)
##########################################################################################################################################

#### II. Assignation des résultats dans les listes correspondantes

for (r in 1:length(results)){  ###de 1 à 9, pour chaque combi de méthode/taxo
#### Network metrics, brown green pathways ####
  clpn[[r]]=results[[r]]$metpoints$connectivity
  clpl[[r]]=results[[r]]$metplots$connectivity
  
  dlpn[[r]]=results[[r]]$metpoints$density
  dlpl[[r]]=results[[r]]$metplots$density
  
  nlpn[[r]]=results[[r]]$metpoints$nb_node
  nlpl[[r]]=results[[r]]$metplots$nb_node
  
  elpn[[r]]=results[[r]]$metpoints$nb_edge
  elpl[[r]]=results[[r]]$metplots$nb_edge
  
  sumbglpn[[r]]=brg_ratios(results[[r]]$brg_point)$sum
  sumbglpl[[r]]=brg_ratios(results[[r]]$brg_plot)$sum
  
  sdbglpn[[r]]=brg_ratios(results[[r]]$brg_point)$sd
  sdbglpl[[r]]=brg_ratios(results[[r]]$brg_plot)$sd
  
  mnbglpn[[r]]=brg_ratios(results[[r]]$brg_point)$mean
  mnbglpl[[r]]=brg_ratios(results[[r]]$brg_plot)$mean
  
  covsbglpn[[r]]=brg_ratios(results[[r]]$brg_point)$nbcovered
  covsbglpl[[r]]=brg_ratios(results[[r]]$brg_plot)$nbcovered
  
  maxsbglpn[[r]]=brg_ratios(results[[r]]$brg_point)$max_dist
  maxsbglpl[[r]]=brg_ratios(results[[r]]$brg_plot)$max_dist

  #### Beta metrics #### 
  bpl[[r]]=results[[r]]$beta_plot$S
  #§bpn[[r]]=results[[r]]$beta_point$S
}

#### III. Visualisation et sauvegarde dans PDF

pdf("Alpha_metrics_0607_sbm_exp_none.pdf")
global_analysis("connectivity","plot",clpl,sh_pl)
global_analysis("connectivity","point",clpn,sh_pn)

#### Density ####
global_analysis("density","plot",dlpl,sh_pl)
global_analysis("density","point",dlpn,sh_pn)

#### Number of nodes and edges ####
global_analysis("Number of vertices","plot",nlpl,sh_pl)
global_analysis("Number of vertices","point",nlpn,sh_pn)

global_analysis("Number of edges","plot",elpl,sh_pl)
global_analysis("Number of edges","point",elpn,sh_pn)

#### Brown-Green Pathways ####
global_analysis("sum of brown green pathways length","point",sumbglpn,sh_pn)
global_analysis("sum of brown green pathways length","plot",sumbglpl,sh_pl)

global_analysis("mean of brown green pathways length","point",mnbglpn,sh_pn)
global_analysis("mean of brown green pathways length","plot",mnbglpl,sh_pl)

global_analysis("sd of brown green pathways length","point",sdbglpn,sh_pn)
global_analysis("sd of brown green pathways length","plot",sdbglpl,sh_pl)

global_analysis("Covered species with brown green pathways","point",covsbglpn,sh_pn)
global_analysis("Covered species with brown green pathways","plot",covsbglpl,sh_pl)

dev.off()

######################################################### Beta ##############################################################################

pdf("Beta_analysis_2007WN_sbm_exp_none.pdf")
#global_analysis("Beta analysis","point",bpn,rep("A",4032))
global_analysis("Beta analysis_WN","plot",bpl,rep("A",360))
dev.off()

sbmsg=lin_mod(bpl[[2]],bpl[[3]])
sbmgf=lin_mod(bpl[[1]],bpl[[2]])

expsg=lin_mod(bpl[[5]],bpl[[6]])
expgf=lin_mod(bpl[[4]],bpl[[5]])

nonesg=lin_mod(bpl[[8]],bpl[[9]])
nonegf=lin_mod(bpl[[7]],bpl[[8]])

beta_wn_rpv=rbind(ldply(sbmsg),ldply(sbmgf),ldply(expsg),ldply(expgf),ldply(nonesg),ldply(nonegf))
combis=rep(c("sbmsg","sbmgf","expsg","expgf","nonesg","nonegf"),each=2)
beta_wn_rpv$ids=combis

######################################################### Analyse de variance ##############################################################################
boxplotting<-function(data,labels,t,xl="plot",yl="metric"){
  df=data.frame(met=data,plots=labels)
  g=boxplot(met~plots,data=df, main=t, xlab=xl, ylab=yl)
  
  print(g)
}

full_combi_metric<-function(l,clabels,rlabels){ ##Valeurs de la métrique pour toutes les combinaisons meth+taxo
  logl=lapply(l, function(x) log10(x)+1)
  mrg=data.frame(do.call("cbind",logl))
  colnames(mrg)=clabels
  mrg$plot=rlabels
  
  return(mrg)
}

collabels=c("SBM_family","SBM_genus","SBM_species","Expert_family","Expert_genus","Expert_species","None_family","None_genus","None_species")
rowlabels=rep(unlist(names(mb_list_t_plot)),each=4)
pdf("anavar_sbmexpnone.pdf")

###Connectivity
c_full=full_combi_metric(clpn,collabels,rowlabels)
boxplotting(c_full$SBM_family,c_full$plot,"Connectivity - SBM - Family")
boxplotting(c_full$SBM_genus,c_full$plot,"Connectivity - SBM - Genus")
boxplotting(c_full$SBM_species,c_full$plot,"Connectivity - SBM - Species")
boxplotting(c_full$Expert_family,c_full$plot,"Connectivity - Expert - Family")
boxplotting(c_full$Expert_genus,c_full$plot,"Connectivity - Expert - Genus")
boxplotting(c_full$Expert_species,c_full$plot,"Connectivity - Expert - Species")
boxplotting(c_full$None_family,c_full$plot,"Connectivity - None - Family")
boxplotting(c_full$None_genus,c_full$plot,"Connectivity - None - Genus")
boxplotting(c_full$None_species,c_full$plot,"Connectivity - None- Species")

###Sum of brown/green pathways length
bg_full=full_combi_metric(sumbglpn,collabels,rowlabels)
boxplotting(bg_full$SBM_family,bg_full$plot,"Brown-green ratios of the sum of pathways length - SBM - Family")
boxplotting(bg_full$SBM_genus,bg_full$plot,"Brown-green ratios of the sum of pathways length - SBM - Genus")
boxplotting(bg_full$SBM_species,bg_full$plot,"Brown-green ratios of the sum of pathways length - SBM - Species")
boxplotting(bg_full$Expert_family,bg_full$plot,"Brown-green ratios of the sum of pathways length - Expert - Family")
boxplotting(bg_full$Expert_genus,bg_full$plot,"Brown-green ratios of the sum of pathways length - Expert - Genus")
boxplotting(bg_full$Expert_species,bg_full$plot,"Brown-green ratios of the sum of pathways length - Expert - Species")
boxplotting(bg_full$None_family,bg_full$plot,"Brown-green ratios of the sum of pathways length - None - Family")
boxplotting(bg_full$None_genus,bg_full$plot,"Brown-green ratios of the sum of pathways length - None - Genus")
boxplotting(bg_full$None_species,bg_full$plot,"Brown-green ratios of the sum of pathways length - None - Species")


dev.off() ##Fermer le fichier pour qu'il soit ouvrable sous Windows

#############################################################################################################
########################### Sortir la liste des espèces qui sautent dans chaque parcelle ####################
ldiffs=list()
cpt=0
for(j in 1:3){
  m=method[j]
  for(i in 1:3){
    tax=taxo[i]
    k=(j-1)*3+i
    netsij=results[[k]]$nets_plot
    for(p1 in 1:16){
      for(p2 in 1:16){
        g1=netsij[[p1]]
        g2=netsij[[p2]]
        
        diff12=toString(setdiff(V(g1)$name,V(g2)$name))
        diff21=toString(setdiff(V(g2)$name,V(g1)$name))
        
        diffs=data.frame(Meth=m,Tax=tax,P1=names(netsij)[p1],P2=names(netsij)[p2],DiffP12=diff12,DiffP21=diff21)
        cpt=cpt+1
        
        ldiffs[[cpt]]=diffs
      }
    }
  }
}

D=data.frame(do.call("rbind",ldiffs))
write.csv2(D,file="tax_compo_diff.csv")
##############################################################################################################

################# LINEAR MODEL #############
dfs=list(length=9)
parcel=rep(names(results[[1]]$nets_plot),each=4)
p=matrix(unlist(strsplit(parcel,split = "\\+|\\-")), ncol=1, byrow=TRUE)
trt=rep(rep(c("0","1"),each=4),4)
sites=c(rep("Mons",16),rep("Lus",24),rep("Theix",24))

for(combi in 1:9){
  n=names(results[[combi]]$nets_point)
  connectivity=clpn[[combi]]
  sumbg=sumbglpn[[combi]]
  df=as.data.frame(cbind(n,connectivity,sumbg))
  parcel=rep(names(results[[1]]$nets_plot),each=4)
  
  df$s=sites
  df$p=p
  df$trt=trt
  
  dfs[[combi]]=df
  
  write.csv2(df,paste(combi,"_lmdata.csv"))
}


# #### Echelle PLOT ####
# yc=c(cssg$R,cesg$R,cnsg$R,csgf$R,cegf$R,cngf$R) ## R²
# x=c(rep(c("SBM","EXPERT","NONE"),2)) 
# z=c(rep("Species_Genus",3), rep("Genus_Family",3))
# 
# alpha_connectivity=alpha_plot(x,yc,z,"Connectivity correlation between taxonomic levels - Plots")
# 
# yd=c(dssg$R,desg$R,dnsg$R,dsgf$R,degf$R,dngf$R) ## R²
# alpha_density=alpha_plot(x,yd,z,titre="Density correlation between taxonomic levels - Plots")
# 
# #### Echelle POINT ####
# p_yc=c(p_cssg$R,p_cesg$R,p_cnsg$R,p_csgf$R,p_cegf$R,p_cngf$R) ## R²
# 
# p_alpha_connectivity=alpha_plot(x,p_yc,z,"Connectivity correlation between taxonomic levels - Sample points")
# 
# p_yd=c(p_dssg$R,p_desg$R,p_dnsg$R,p_dsgf$R,p_degf$R,p_dngf$R) ## R²
# p_alpha_density=alpha_plot(x,p_yd,z,titre="Density correlation between taxonomic levels - Sample points")
# 
# 

# ### Echelle plot ###
# sx=c(results[[3]]$metplots$connectivity,results[[6]]$metplots$connectivity,results[[9]]$metplots$connectivity)
# sz=rep(c("SBM","EXPERT","NONE"),each=length(results[[3]]$metplots$connectivity))
# sy=c(results[[2]]$metplots$connectivity,results[[5]]$metplots$connectivity,results[[8]]$metplots$connectivity)
# 
# alpha_smooth(sx,sy,sz,"Connectivity")
# 
# dsx=c(results[[3]]$metplots$density,results[[6]]$metplots$density,results[[9]]$metplots$density)
# dsz=rep(c("SBM","EXPERT","NONE"),each=length(results[[3]]$metplots$density))
# dsy=c(results[[2]]$metplots$density,results[[5]]$metplots$density,results[[8]]$metplots$density)
# 
# alpha_smooth(dsx,dsy,dsz,"Density")
# 
# ### Echelle point ###
# shapes=rep(c(rep("Mons",16),rep("Lusignan",24),rep("Theix",24)),3)
# p_sx=c(results[[3]]$metpoints$connectivity,results[[6]]$metpoints$connectivity,results[[9]]$metpoints$connectivity)
# p_sz=rep(c("SBM","EXPERT","NONE"),each=length(results[[3]]$metpoints$connectivity))
# p_sy=c(results[[2]]$metpoints$connectivity,results[[5]]$metpoints$connectivity,results[[8]]$metpoints$connectivity)
# 
# alpha_smooth(p_sx,p_sy,p_sz,"Connectivity",shapes)+
#   coord_trans(x = "log10", y="log10")
# 
# p_dsx=c(results[[3]]$metpoints$density,results[[6]]$metpoints$density,results[[9]]$metpoints$density)
# p_dsz=rep(c("SBM","EXPERT","NONE"),each=length(results[[3]]$metpoints$density))
# p_dsy=c(results[[2]]$metpoints$density,results[[5]]$metpoints$density,results[[8]]$metpoints$density)
# 
# alpha_smooth(p_dsx,p_dsy,p_dsz,"Density",shapes)+
#   coord_trans(x = "log10", y="log10")
# ########################################
# 
# ### SBM
# ### Echelle plot
# beta_plot_fam=results[[1]]$beta_plot$S
# beta_plot_gen=results[[2]]$beta_plot$S
# beta_plot_spec=results[[3]]$beta_plot$S
# 
# plot_sg=compare_beta(beta_plot_spec,beta_plot_gen,"Species","Genus","Plot beta diversity")
# plot_gf=compare_beta(beta_plot_gen,beta_plot_fam,"Genus","Family","Plot beta diversity")
# 
# ### Echelle point
# beta_point_fam=results[[1]]$beta_point$S
# beta_point_gen=results[[2]]$beta_point$S
# beta_point_spec=results[[3]]$beta_point$S
# 
# point_sg=compare_beta(beta_point_spec,beta_point_gen,"Species","Genus","Point beta diversity")
# point_gf=compare_beta(beta_point_gen,beta_point_fam,"Genus","Family","Point beta diversity")
# 
# ######################################################################################################################################## 
# mybeta<-results[[1]]$beta_point
# hist(mybeta$S)
# mybeta$S1<- cut(mybeta$S, breaks = c(0.75,0.82, 0.9, 1),right = FALSE)
# ggplot(data = mybeta, aes(x=i, y=j)) + 
#   geom_tile(aes( fill = S), colour = "white")#+
# #scale_fill_brewer(palette = "PRGn")
# 
# #network_betaplot(ntwk_list_plot[[1]], ntwk_list_plot[[12]]) ### Vert premi?re parcelle (r?seau) exclusivement, bleu 2?me parcelle, gris les deux
# 


##### Unit tests #####

empties=list()
for(x in results){
  for(i in 1:64){
    if(x$metpoints$nb_node[i]==0){
      empties=append(empties,paste(x$meth,x$tax,i,sep=" | "))
    }
  }
}
