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

###4. Build graph
taxonomy<-as.Node(names_tree)
#kg<-as.igraph(taxonomy,vertexAttributes = c("key","occurrenceId","verbatimScientificName"))

###5. Complete  trophic lists with info inherited from parents
Complete_troph<-function(node){  
  ##Traverses the taxonomic tree and on each node, appends to trophic columns the values of all their ancestors
  node$res_full=array(unlist(node$Get('resources', traversal = "ancestor")))
  node$cons_full=node$res_full[!is.na(node$res_full)]
  node$cons_full=array(unlist(node$Get('consumers', traversal = "ancestor")))
  node$cons_full=node$cons_full[!is.na(node$cons_full)]
  node$host_full=array(unlist(node$Get('hosts', traversal = "ancestor")))
  node$host_full=node$host_full[!is.na(node$host_full)]
  node$paras_full=array(unlist(node$Get('parasites', traversal = "ancestor")))
  node$paras_full=node$paras_full[!is.na(node$paras_full)]
  
  return(node)
}

taxonomy$Do(Complete_troph, filterFun = isNotRoot)
print(taxonomy,"key","res_full","host_full","paras_full","cons_full",limit=10)


###6. Build interaction network
###Traverse the tree and create an igraph node with verbatimSN and add links
nodes=Traverse(taxonomy,filterFun = Relev_node)

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

######### Create vertices ############
InitNet<-function(nodes){
  ###Transform to igraph node
  int_net<-make_empty_graph()
  for(node in nodes){
    v=vertex(key=node$key,name=as.character(node$verbatimScientificName))
    int_net=int_net+v
  }
  return(int_net)
}

network=InitNet(nodes)

######### Create edges ############
AddLinks<-function(nodes,net){
  for(node in nodes){
    res_keys=node$res_full
    nv=V(net)[key==node$key]
    for(neigh in res_keys){
      if(!is.na(neigh)){
        print(neigh)
        neighv=V(net)[key==neigh]
        if(length(neighv)==0){
          neighv=vertex(key=neigh,name="?")
          net=net+neighv
        }
        e=edge(nv,neighv)
        net=net+e
      }
    }
  }
}

full_network=AddLinks(nodes,network)
