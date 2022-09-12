get_simplified_topology<-function (subnodes, network,links){
  
  # browser()
  #source("./data/functions/F06_get_link.R")
  
  
  if(all(subnodes!='')){
    subnodes<-unlist(subnodes)
    
    network<-igraph::simplify(network, remove.multiple = T, remove.loops = F)
    
    num_nodes <-length(subnodes)
    p <- 1:num_nodes
    #NL_positions <- num_nodes*(p-1)+p 
    positions <- 1:num_nodes^2
    #positions <- positions[!is.element(positions,NL_positions)]
    
    combi_nodes <- gtools::permutations(num_nodes, 2, 1:num_nodes, repeats.allowed=T)
    combi_nodes <- combi_nodes[ ,c(2,1)]
    

    tuple <-rep(0,num_nodes^2)
    for (i in 1:nrow(combi_nodes)){
      tmp_nodes <- p[-combi_nodes[i,]] # 
      tmp_names <- subnodes[tmp_nodes]
      tmp_net <- igraph::delete_vertices(network, tmp_names)
      tmp<-get_link(tmp_net,subnodes[combi_nodes[i,1]],subnodes[combi_nodes[i,2]],links, full= F)
      if(!is.na(tmp)){
        
        tuple[i] <- tmp}
    }
    
    
    return(tuple)}
  
}
