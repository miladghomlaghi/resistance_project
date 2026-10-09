

# library(igraph)
library(prodlim) #for row match
# #library(DiagrammeR)
# #library(DiagrammeRsvg)
# library(foreach)
# library(doParallel)
# library(readxl)
# library(tidyverse)
# library(tictoc)

source("./R/Functions/F02_get_sign.R")
source("./R/Functions/F03_get_link.R")

start_node <- 1
end_node <- 25
data_to_match <- 'signor_PathwaysPlus'
topos_to_match_name <- 'robust_minimal_4Nodes_Seq'
num_nodes <- 4
nodesAreTargets <- c(1,2)


load(paste0("./Data/matching/",data_to_match,".RData", collapse = ''))
load(paste0("./Data/matching/",num_nodes,"Nodes/",topos_to_match_name, ".RData", collapse = ''))
load(paste0("./Data/matching/","Targets.RData", collapse = ''))
load(paste0("./Data/topologies/",num_nodes,"Nodes/topos_all.RData", collapse = ''))

#topos_to_match<-topos8rob_to_match
#save(topos_to_match,file='topos31_3Nodes.RData')


targets <- nodes[is.element(nodes,targets)]
targets_num<-match(targets, nodes)
#save(topos_all, file="topologies_4Nodes_binary")
links$ENTITYA<-gsub("/",";",links$ENTITYA) 
links$ENTITYB<-gsub("/",";",links$ENTITYB)
links$ENTITYA<-gsub(",",":",links$ENTITYA) 
links$ENTITYB<-gsub(",",":",links$ENTITYB)
nodes<-gsub("/",";",nodes)
nodes<-gsub(",",":",nodes)
targets <- gsub(",",":",targets)
targets<-gsub("/",";",targets)
vertex_attr(net)$name<-gsub("/",";",vertex_attr(net)$name)
vertex_attr(net)$name<-gsub(",",":",vertex_attr(net)$name)

pracma::tic()
final_results<-c(0,"none")
#final_resultsNeg<-c(0,"none")

#num_nodes <-sqrt(ncol(topos_all))
p <- 1:num_nodes
NL_positions <- num_nodes*(p-1)+p 
positions <- 1:num_nodes^2
positions <- positions[!is.element(positions,NL_positions)]
combi_nodes <- gtools::permutations(num_nodes, 2, 1:num_nodes, repeats.allowed=F)
combi_nodes <- combi_nodes[order(combi_nodes[,2]),]


for (v in start_node:end_node){
  
  num <- which(nodes==targets[v])
  range <- 1:length(nodes)
  
  numCores<-bigstatsr::nb_cores()
  doParallel::registerDoParallel(numCores)
  
  all_combi<- gtools::permutations(length(nodes)-1, num_nodes-1, range[-num], repeats.allowed=F) 
  if (length(nodesAreTargets)>1){
  all_combi<-all_combi[is.element(all_combi[,1],targets_num),]
  }

  
  #numCores<-parallel::detectCores(all.tests = FALSE, logical = TRUE)
  
  results <- foreach::foreach(c=1:nrow(all_combi),.combine = "rbind") %dopar% {
  god_com <- ""
 # for (c in 1:100)  {
   #(igraph)
    #library(prodlim)
    
    #god_comNeg <- 0
    combi<- replicate(num_nodes,0)
    if (length(nodesAreTargets)>1){
      combi[nodesAreTargets[1]] <- num
      combi[nodesAreTargets[2]] <- all_combi[c,1]
      combi[p[!is.element(p,nodesAreTargets)]] <- all_combi[c,2:(num_nodes-1)] 
    }else{
      combi[nodesAreTargets[1]] <- num
      combi[p[!is.element(p,nodesAreTargets)]] <- all_combi[c,3:(num_nodes-1)]
    }
    
    
    
    
    
    combi <- matrix(combi,nrow=1,ncol = num_nodes)
    combi_names<-nodes[combi]
    topo <-rep(0,num_nodes^2)
    for (i in 1:nrow(combi_nodes)){
      tmp_nodes <- p[-combi_nodes[i,]] # 
      tmp_names <- combi_names[tmp_nodes]
    tmp_net <- igraph::delete_vertices(net, tmp_names)
      topo[positions[i]] <- get_link(tmp_net,combi_names[combi_nodes[i,1]],combi_names[combi_nodes[i,2]],links)
    }
    if (!is.na(prodlim::row.match(topo,topos_to_match))){
      god_com<-c(god_com,paste0(combi_names,collapse = ","))
    }
    # if (!is.na(row.match(topoNeg,topos8rob_to_match))){
    #   god_comNeg<-c(god_comNeg,c)
   
    data.frame(god_com)
    }
   
  
    #}
  #results<-
 
  
   stopImplicitCluster()
  
  
  if (length(which(results!=""))==0){
    
    final_results<-rbind(final_results,c(num,"none"))
  }else{
    final_results<-rbind(final_results,c(num,paste(results[which(results!=""),],collapse="/")))
    
  }
}


final_results<-final_results[-1,]
#final_resultsNeg<-final_resultsNeg[-1,]
save(object = final_results, file = paste0("run",start_node,"to",end_node,"_",topos_to_match_name,"_nodes",num_nodes,".RData"))
# save(object = results, file = paste0("run",start_node,"to",end_node,"_",data_to_match,"_nodes",num_nodes,".RData"))

#save(object = final_resultsNeg, file = "run1to20Neg")
pracma::toc()

#```
