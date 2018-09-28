library(mixer)
library(ggplot2)

n=50
meth="variational"

mixs=mixbayes
mixg=mixbayesg
mixf=mixbayesf

setwd(".")
dataf=read.csv2("adjacence_esp_bin_family.csv",header=TRUE,row.names=1)

mixbayesf=mixer(dataf,qmin=1,qmax=n,method = meth, nbiter=50, fpnbiter=25, improve=TRUE)

mixbayes=mixbayesf
#### ICL ####
icl=vector(mode='numeric',length=n)
index=1:n
write.csv2(mixbayes$nnames[[1]],file="names_family.csv")
for (i in index){
  ### get ICL per number of groups
  icl[i]=mixbayes$output[i][1][[1]]$criterion
  ### get connectivity matrix per number of groups
  m=data.matrix(mixbayes$output[i][1][[1]]$Pis)
  write.csv2(m,file=paste(c("OUT_R_MIXNET_FAMILY/Q",i,"groups.csv"),collapse="_"))
  ### get partitioning per number of groups
  taus=as.matrix(mixbayes$output[i][1][[1]]$Taus)
  #p=t(as.matrix(1:i))%*%taus
  #out=t(p)
  #out$spec=mixbayes$nnames[[1]]
  out=data.frame(row.names = mixbayes$nnames[[1]])
  #out$species=mixbayes$nnames[[1]]
  lv=split(taus, rep(1:ncol(taus), each = nrow(taus)))
  out$labels=unlist(lapply(lv,FUN = function(x) which.max(x)))
  write.csv2(out,file=paste(c("OUT_R_MIXNET_FAMILY/Q",i,"labels.csv"),collapse="_"))
}
plot(icl~index)
qoptim=which.max(icl)
maxicl=max(icl)
write.csv2(as.data.frame(icl),file="family_icl.csv")
