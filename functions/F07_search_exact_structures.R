
find_exact_structures <- function(topology, user_node, user_position, nodes, links, network) {
  
  #  topology<-c(-1,1,1,0,0,1,0,0,0)
  #   user_node<-"EGFR"
  #    user_position <-1
  #   nodes<-nodes[[3]]
  #    links<-links[[3]] 
  #   network<-nets[[3]]
  
  dist<-5
  positions <- 1:length(topology)
  num_nodes <- sqrt(length(topology))
  p <- 1:num_nodes
  nodes<-sort(nodes)
  #NL_positions <- num_nodes*(p-1)+p
  
  #topology[NL_positions]<-0
  
  
  #user_node_num <- which(stringi::stri_cmp_eq(user_node,nodes))
  user_node_num <- which(is.element(nodes,as.character(user_node)))
  
  #final_resultsNeg<-c(0,"none")
  
  #num_nodes <-sqrt(ncol(topos_all))
  
  
  #positions <- positions[!is.element(positions,NL_positions)]
  combi_nodes <- gtools::permutations(num_nodes, 2, 1:num_nodes, repeats.allowed=T)
  combi_nodes <- combi_nodes[ ,c(2,1)]
  
  
  if(length(nodes)>300){
    near_net<-igraph::ego(network,dist, nodes=as.character(user_node), mode = "out", mindist = 0)
    network<-igraph::induced_subgraph(network,unlist(near_net))
    
    nodes<-igraph::as_ids(igraph::V(network))
    nodes<-sort(nodes)
  }
  
  
  range <- 1:length(nodes)
  
  
  numCores<-bigstatsr::nb_cores()
  doParallel::registerDoParallel(cl <- makeCluster(numCores))
  
  all_combi<- gtools::permutations(length(nodes)-length(user_node_num), num_nodes-1, range[-user_node_num], repeats.allowed=F)
  #numCores<-parallel::detectCores(all.tests = FALSE, logical = TRUE)
  p1 <- p[-user_position][1]
  p2 <- p[-user_position][2]
  
  results <- foreach::foreach(r=1:nrow(all_combi),.combine = "rbind") %dopar% {
    source("./data/functions/F02_get_link.R")
    
    #    results <-matrix(rep("none",num_nodes),1,num_nodes)
    
    #   for (r in 1:nrow(all_combi)){
    
    combi <- rep(0, num_nodes)
    combi[user_position] <- user_node_num
    combi[p1] <- all_combi[r, ][1]
    combi[p2] <- all_combi[r, ][2]
    combi <- matrix(combi, nrow = 1, ncol = num_nodes)
    combi_names <- nodes[combi]
    topo <- rep(0, num_nodes ^ 2)
    
    for (i in 1:nrow(combi_nodes)) {
      tmp_nodes <- p[-combi_nodes[i, ]] #
      tmp_names <- combi_names[tmp_nodes]
      tmp_net <- igraph::delete_vertices(network, tmp_names)
      tmp_link <-
        get_link(tmp_net, combi_names[combi_nodes[i, 1]], combi_names[combi_nodes[i, 2]], links, full =
                   F)
      if (!is.na(tmp_link)) {
        topo[positions[i]] <- tmp_link
      }
    }
    if (sum(topo == topology) == num_nodes ^ 2) {
      combi_names
    }
  }
  
  stopCluster(cl)
  results<-as.data.frame(results)
  
  if(nrow(results)<1){
    return(rbind(results,matrix(rep("none",num_nodes),1,num_nodes)))
  } else{
    return(results)
  }
  
  
}
