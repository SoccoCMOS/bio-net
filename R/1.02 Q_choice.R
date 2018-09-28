library(mixer)
library(ggplot2)

n=50
eps=0.05  ##Tolérance en dessous du nombre de groupe optimal

taxof=c("family","genus","species")
repos=c("../MIXNET_TAXO/OUT_R_MIXNET_FAMILY","../MIXNET_TAXO/OUT_R_MIXNET_GENUS","../MIXNET_TAXO/OUT_R_MIXNET_SPEC")

iclfiles=c("family_icl.csv","genus_icl.csv","spec_icl.csv")

outpu_dir="../data/input/taxo_metabar/adj/sbm"

setwd(".")

qopt=vector(mode='numeric',length=3)
for (taxo in 1:3){
  rep=repos[taxo]
  fil=iclfiles[taxo]
  
  icl_file_path=paste(rep,fil,sep="/")
  
  icltab=data.matrix(read.csv2(icl_file_path,header=TRUE,row.names=1))
  
  iclmax=max(icltab)
  qmax=which.max(icltab)
  i=qmax
  
  target=(1+eps)*iclmax
  
  while(icltab[i]>=target & i<=n){  ##> car c'est négatif
    i=i+1
  }
  qopt[taxo]=i
  
  adj_mat_file=paste(rep,"/","Q_",i,"_groups.csv",sep="")
  adj_mat_out=paste(outpu_dir,"/",taxof[taxo],".csv",sep="")
  
  adj_mat=read.csv2(adj_mat_file,header=TRUE,row.names=1)
  row.names(adj_mat)=colnames(adj_mat)
  write.csv2(adj_mat,adj_mat_out)
}



