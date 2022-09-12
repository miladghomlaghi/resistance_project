get_size_logic<-function (proteins, net,links, size){
  
  source("./data/functions/F05_get_sign.R")
  source("./data/functions/F06_get_link.R")
  
  proteins<-unlist(proteins)
  
  
  
  num_nodes <-length(proteins)
  p <- 1:num_nodes
  NL_positions <- num_nodes*(p-1)+p 
  positions <- 1:num_nodes^2
  positions <- positions[!is.element(positions,NL_positions)]
  combi_nodes <- gtools::permutations(num_nodes, 2, 1:num_nodes, repeats.allowed=F)
  combi_nodes <- combi_nodes[order(combi_nodes[,2]),]
  
  
  
  
  topo <-rep(T,num_nodes^2)
  for (i in 1:nrow(combi_nodes)){
    tmp_nodes <- p[-combi_nodes[i,]] # 
    tmp_names <- proteins[tmp_nodes]
    tmp_net <- igraph::delete_vertices(net, tmp_names)
    topo[positions[i]] <- get_size(tmp_net,proteins[combi_nodes[i,1]],proteins[combi_nodes[i,2]],links,size)
  }
  

  
  return( sum(topo) == num_nodes^2 )
  
}
