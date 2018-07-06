# -----------------------------------------------------------------------------
# PROJECT:
#    BISE : Linking soil foodweb and soil functions dynamics under
#    agricultural practices (des)intensification and different input levels
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# STEP:
#    0.1   Generating species x species adjacency matrices (binary, weighted by proba & count)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# NOTES:
#  
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# LIBRARIES:
library(igraph)
library(RNewsflow)
library(reshape2)
# -----------------------------------------------------------------------------

##### Params ####

#################

setwd("../data/input/taxo_metabar/adj")

# -----------------------------------------------------------------------------
# data load
# -----------------------------------------------------------------------------
#Adjacency matrix
adj     <- read.csv("adjacence_esp_bin_spec.csv", sep = ";", h = T, row.names=1) 
adj     <- as.matrix(adj)
spec_names <- rownames(adj)
colnames(adj) <- spec_names
#FW Codes
codesFW <- read.csv("../species_species_metabar_raw.csv", sep = ";", h = T) 
u_codes <- unique(codesFW[,c("retained_tax","codeFW")])
codes <- subset(u_codes, u_codes$retained_tax %in% spec_names)
# -----------------------------------------------------------------------------


longadj=melt(as.matrix(adj), na.rm = TRUE)
el1=merge(x=longadj,y=codes,by.x="Var2",by.y = "retained_tax",all=TRUE,sort=TRUE)
el2=merge(x=el1,y=codes,by.x="Var1",by.y = "retained_tax",suffixes=c("_2","_1"),sort=TRUE)

#Aggregation par codeFW
agg_adj_list=aggregate(el2$value,by=list(el2$codeFW_1,el2$codeFW_2),FUN=function(x) ifelse(sum(x)>0,1,0))

#Subset des interactions égales à 1
#agg_adj_list_filt=subset(agg_adj_list,agg_adj_list$x>0)

agg_adj_mat=dcast(agg_adj_list,Group.1~Group.2,value.var="x",fill=0)
rownames(agg_adj_mat) <- agg_adj_mat[,1]
agg_adj_mat<-agg_adj_mat[,-1]

write.csv2(agg_adj_mat,"species_expert_agg.csv")

