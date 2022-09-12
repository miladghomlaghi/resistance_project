get_subnetwork <- function (nodes, network) {
  warning_state <- getOption("warn")
  options(warn = -1)
  
  # if (!is.vector(subnodes)){
  #   
  #     subnodes<- c(as.character(subnodes[1,1]),
  #              as.character(subnodes[1,2]),
  #              as.character(subnodes[1,3]))
  # 
  # }
  
  nodes <- sort(as.character(nodes))
  num_nodes <-length(nodes)
  p <- 1:num_nodes
  #NL_positions <- num_nodes*(p-1)+p 
  #positions <- 1:num_nodes^2
  #positions <- positions[!is.element(positions,NL_positions)]
  combi_nodes <- gtools::permutations(num_nodes, 2, 1:num_nodes, repeats.allowed=T)
  combi_nodes <- combi_nodes[order(combi_nodes[,2]),]
  
   #topo <-rep(0,num_nodes^2)
  # for (i in 1:nrow(combi_nodes)){
  #   tmp_nodes <- p[-combi_nodes[i,]] # 
  #   tmp_names <- proteins[tmp_nodes]
  #   tmp_net <- igraph::delete_vertices(network, tmp_names)
  #   topo[positions[i]] <- get_link(tmp_net,proteins[combi_nodes[i,1]],proteins[combi_nodes[i,2]],links)
  # }

  tmp_nodesName <-""
  for (i in 1:nrow(combi_nodes)){
    tmp_nodes <- p[-combi_nodes[i,]] # 
    tmp_names <- nodes[tmp_nodes]
    #tmp_net <- igraph::delete_vertices(network, tmp_names)
    #topo[positions[i]] <- get_link(tmp_net,proteins[combi_nodes[i,1]],proteins[combi_nodes[i,2]],links)
    L_net <- igraph::delete_vertices(network, tmp_names)
    L_short <- names(unlist(igraph::shortest_paths(L_net,nodes[combi_nodes[i,1]],nodes[combi_nodes[i,2]])[[1]]))
    #L_shortB<-names(unlist(igraph::shortest_paths(L_net,subnodes[2],subnodes[1])[[1]]))
    
    tmp_nodesName <- c(tmp_nodesName,L_short)
  }
  tmp_nodesName<-tmp_nodesName[-1]
# L1_net <- igraph::delete_vertices(network, subnodes[3])
# L1_short<-names(unlist(igraph::shortest_paths(L1_net,subnodes[1],subnodes[2])[[1]]))
# L1_shortB<-names(unlist(igraph::shortest_paths(L1_net,subnodes[2],subnodes[1])[[1]]))
# proteins[combi_nodes[i,2]]
# 
# L2_net <- igraph::delete_vertices(network, subnodes[2])
# L2_short<-names(unlist(igraph::shortest_paths(L2_net,subnodes[1],subnodes[3])[[1]]))
# L2_shortB<-names(unlist(igraph::shortest_paths(L2_net,subnodes[3],subnodes[1])[[1]]))
# 
# 
# 
# L3_net <- delete_vertices(network, subnodes[1])
# L3_short<-names(unlist(igraph::shortest_paths(L3_net,subnodes[2],subnodes[3])[[1]]))
# L3_shortB<-names(unlist(igraph::shortest_paths(L3_net,subnodes[3],subnodes[2])[[1]]))


#sub_nodes<-unique(c(L1_short,L1_shortB,L2_short,L2_shortB,L3_short,L3_shortB))
sub_nodes<-unique(tmp_nodesName)
sub_V1_Net<-igraph::induced_subgraph(network,sub_nodes)

options(warn = warning_state)
return(sub_V1_Net)

}
