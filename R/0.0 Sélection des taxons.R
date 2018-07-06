# -----------------------------------------------------------------------------
# PROJECT:
#    BISE : Linking soil foodweb and soil functions dynamics under
#    agricultural practices (des)intensification and different input levels
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# STEP:
#    0.0   Selection of taxa at parameters level
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# NOTES:
#  
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# LIBRARIES:
library(pracma)
# -----------------------------------------------------------------------------

##### Params ####
input="metabar_raw.csv"

### Possible values of high and low: phylum, class, order, family, genus, species ###
tax_levels=c("phylum","class","order", "family", "genus", "species")

high="family"
low="family"

#################

setwd("../data/input")

mb_raw=read.csv2(input,sep=";",h=T)

### Control integrity of input parameters ###
highctr=paste(high,"name",sep="_")
lowctr=paste(low,"name",sep="_")

high_ind=which(tax_levels==high)
low_ind=which(tax_levels==low)

if(high_ind<low_ind){
  tmp=high_ind
  high_ind=tmp
  low_ind=tmp
}

if (!((highctr %in% colnames(mb_raw))&(highctr %in% colnames(mb_raw)))){
  print("Taxonomic level not provided in dataset ")
}else {
  ### Get target columns
  target_col_names=lapply(tax_levels[high_ind:low_ind],function(x) paste(x,"name",sep="_"))
  target_cols=mb_raw[unlist(target_col_names)]

  ### Construction of retained taxa ###
  retained<-vector(mode="character",length=nrow(target_cols))
  
  todrop<-list()
  
  for (i in 1:nrow(target_cols)){
    notna=which(!is.na(target_cols[i,]))
    if (length(notna)>0) {
      retained[i] <- as.character(target_cols[i,tail(notna,n=1)])
    }
    else{
      retained[i] <-NA
      todrop <- append(todrop,i)
    }
  }
  
  ### Filtered taxa 
  #removed<-mb_raw[unlist(todrop),]
  mb<-mb_raw[-unlist(todrop),]
  ret<-retained[-unlist(todrop)]
  
  ### Add to original matrix
  mb$retained_tax=ret
  
  tax_conf <- mb[,c("retained_tax","codeFW","µhab_surf","µhab_subsurf","µhab_soil")]
  
  conflict_codeFW=list()
  conflict_hab=list()
  
  ### Electing representative codeFW => aggregation => TODO: Conflicts: keep all vs keep mode
  for (t in unique(ret)){
    indices=which(tax_conf$retained_tax %in% c(t))
    
    ### Aggregate codeFW
    troph_codes=tax_conf[indices,]$codeFW ## Conflicting trophic codes.
    if(length(unique(troph_codes))>1){
      conflict_codeFW=append(conflict_codeFW,t)
    }
    mod_value=Mode(troph_codes)
    tax_conf[indices,]$codeFW<-rep(mod_value,length(indices))
    
    ### Aggregate µhabitat => µhab_surf is = 1 if any (at least one) row with retained_taxa=t has µhab_surf=1 (same for subsurf, soil)
    hab=tax_conf[indices,c("retained_tax","µhab_surf","µhab_subsurf","µhab_soil")]
    agg=aggregate(hab[,c("µhab_surf","µhab_subsurf","µhab_soil")],by=list(hab$retained_tax),function(x) ifelse(sum(x)>0,1,0)) ##Logical OR
    
    tax_conf[indices,]$µhab_surf<-rep(agg$µhab_surf,length(indices))
    tax_conf[indices,]$µhab_subsurf<-rep(agg$µhab_subsurf,length(indices))
    tax_conf[indices,]$µhab_soil<-rep(agg$µhab_soil,length(indices))
  }
  
  mb[,c("retained_tax","codeFW","µhab_surf","µhab_subsurf","µhab_soil")]<-tax_conf
  
  
  ### Save to output file
  file_out=paste(high,low,input,sep="_")  
  
  write.csv2(mb,file=file_out)
}


