get_simplified_topology<-function (proteins, net,links){

  
  source("./R/Functions/F02_get_sign.R")
  source("./R/Functions/F03_get_link.R")
  
  # if (!is.vector(proteins)){
  #   
  #   proteins<- c(as.character(proteins[1,1]),
  #                as.character(proteins[1,2]),
  #                as.character(proteins[1,3]))
  #   
  # }
  
  num_nodes <-length(proteins)
  p <- 1:num_nodes
  NL_positions <- num_nodes*(p-1)+p 
  positions <- 1:num_nodes^2
  positions <- positions[!is.element(positions,NL_positions)]
  combi_nodes <- gtools::permutations(num_nodes, 2, 1:num_nodes, repeats.allowed=F)
  combi_nodes <- combi_nodes[order(combi_nodes[,2]),]
  
  
  
  
  topo <-rep(0,num_nodes^2)
  for (i in 1:nrow(combi_nodes)){
    tmp_nodes <- p[-combi_nodes[i,]] # 
    tmp_names <- proteins[tmp_nodes]
    tmp_net <- igraph::delete_vertices(net, tmp_names)
    topo[positions[i]] <- get_link(tmp_net,proteins[combi_nodes[i,1]],proteins[combi_nodes[i,2]],links)
  }
  
  # topo <-rep(0,9)
  # tmp_net <- igraph::delete_vertices(net, proteins[3])
  # topo[2] <- get_link(tmp_net,proteins[2],proteins[1],links)
  # tmp_net <- igraph::delete_vertices(net, proteins[2])
  # topo[3] <- get_link(tmp_net,proteins[3],proteins[1],links)
  # tmp_net <- igraph::delete_vertices(net, proteins[3])
  # topo[4] <- get_link(tmp_net,proteins[1],proteins[2],links)
  # tmp_net <- igraph::delete_vertices(net, proteins[1])
  # topo[6] <- get_link(tmp_net,proteins[3],proteins[2],links)
  # tmp_net <- igraph::delete_vertices(net, proteins[2])
  # topo[7] <- get_link(tmp_net,proteins[1],proteins[3],links)
  # tmp_net <- igraph::delete_vertices(net, proteins[1])
  # topo[8] <- get_link(tmp_net,proteins[2],proteins[3],links)
  
return(topo)
}