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

# -----------------------------------------------------------------------------

##### Params ####
input="metabar_raw.csv"

### Possible values of high and low: phylum, class, order, family, genus, species ###
tax_levels=c("phylum","class","order", "family", "genus", "species")

high="family"
low="genus"

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
  mb$ret_taxa=ret
  
  ### Check for conflicting codeFW
  tocheck <- mb[,c("ret_taxa","codeFW")]
  
  ### Save to output file
  file_out=paste(high,low,input,sep="_")  
  
  write.csv2(mb_raw,file=file_out)
}


