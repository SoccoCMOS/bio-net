library(igraph)
library(data.tree)

########### Knowledge graph ###################

###1. Taxonomy 
names<-read.csv2("0.gbif_normalized_names.csv")
names$occurrenceId<-1:dim(names)[1]

taxo_levels<-c("kingdom","phylum","class","order","family","genus","species")

names_tree<-names[append(taxo_levels,c("key","occurrenceId","verbatimScientificName"))]
names_tree$pathString <- paste("root2", 
                               names$kingdom, 
                               names$phylum,
                               names$class,
                               names$order,
                               names$family,
                               names$genus,
                               names$species,
                               sep = "/")

###2. Trophic properties
troph<-read.csv2("1.interactions_list.csv")

###3. Filling strict consumer and resource list
##Assign to each taxa list of consumers and parasites + list of resources and hosts
names_tree$resources<-lapply(names_tree$key,get_out,"eats")
names_tree$hosts<-lapply(names_tree$key,get_out,"parasiteOf")
names_tree$consumers<-lapply(names_tree$key,get_in,"eats")
names_tree$parasites<-lapply(names_tree$key,get_in,"parasiteOf")

get_out<-function(k,type="eats"){
  res=subset(troph,consumer_gbif_key==k & interaction_type==type)$resource_gbif_key
  return(res[!is.na(res)])
}

get_in<-function(k,type="eats"){
  res=subset(troph,resource_gbif_key==k & interaction_type==type)$consumer_gbif_key
  return(res[!is.na(res)])
}

###4. Build tree
taxonomy<-as.Node(names_tree)

###5. Complete  trophic lists with info inherited from parents
Complete_troph<-function(node){  
  ##Traverses the taxonomic tree and on each node, appends to trophic columns the values of all their ancestors
  node$res_full=array(unlist(node$Get('resources', traversal = "ancestor")))
  node$res_full=node$res_full[!is.na(node$res_full)]
  node$cons_full=array(unlist(node$Get('consumers', traversal = "ancestor")))
  node$cons_full=node$cons_full[!is.na(node$cons_full)]
  node$host_full=array(unlist(node$Get('hosts', traversal = "ancestor")))
  node$host_full=node$host_full[!is.na(node$host_full)]
  node$paras_full=array(unlist(node$Get('parasites', traversal = "ancestor")))
  node$paras_full=node$paras_full[!is.na(node$paras_full)]
  
  return(node)
}

taxonomy$Do(Complete_troph, filterFun = isNotRoot)
print(taxonomy,"key","host_full","paras_full","res_full","cons_full",limit=20)

###6. Generate edge list
##Traversal 
Relev_node<-function(x){
  if(isNotRoot(x)){
    k=x$key
    if(is.null(k)){
      return(FALSE)
    }else{
      if(is.na(k)){
        return(FALSE)
      }else{
        return(TRUE)
      }
    }
  }else{
    return(FALSE)
  }
}
nodes=Traverse(taxonomy,filterFun = Relev_node)

Generate_Edges<-function(nodes){
  df <- data.frame(matrix(ncol = 3, nrow = 0))
  x <- c("consumer", "resource", "type")
  colnames(df) <- x
  
  for(node in nodes){
    res=node$res_full
    cons=rep(node$key,length(res))
    type=rep("eats",length(res))
    df=rbind(df,data.frame("consumer"=cons,"resource"=res,"type"=type))
    
    cons=node$cons_full
    res=rep(node$key,length(cons))
    type=rep("eats",length(cons))
    df=rbind(df,data.frame("consumer"=cons,"resource"=res,"type"=type))
    
    host=node$host_full
    paras=rep(node$key,length(host))
    type=rep("parasiteOf",length(host))
    df=rbind(df,data.frame("consumer"=paras,"resource"=host,"type"=type)) 
    
    paras=node$paras_full
    host=rep(node$key,length(paras))
    type=rep("parasiteOf",length(paras))
    df=rbind(df,data.frame("consumer"=paras,"resource"=host,"type"=type))     
  }
  
  return(df)
}

edf=Generate_Edges(nodes)

###Add names to taxa
nm=names[c("key","verbatimScientificName")]
nm_cons=troph[c("consumer_gbif_key","consumer_name")]
colnames(nm_cons)=c("key","verbatimScientificName")
nm_res=troph[c("resource_gbif_key","resource_name")]
colnames(nm_res)=c("key","verbatimScientificName")
nm_troph=rbind(nm_cons,nm_res)

###Check for conflictual associations
key_names=unique(rbind(nm,nm_troph))

Get_Name<-function(k){
  opt=subset(key_names,key==k)$verbatimScientificName
  # if(length(opt)>1){
  #   print(k)
  #   print("Warning possible conflicting keys")
  # }
  name=as.character(opt[1])
  return(name)
}

edf$consumer_name=lapply(edf$consumer,Get_Name)
edf$resource_name=lapply(edf$resource,Get_Name)
