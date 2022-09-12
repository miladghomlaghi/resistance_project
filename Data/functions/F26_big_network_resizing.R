big_network_resizing <- function(network,user_node){
near_net <-
  igraph::ego(
    network,
    dist,
    network_proteins = as.character(user_node),
    mode = "out",
    mindist = 0
  )
network <- igraph::induced_subgraph(network, unlist(near_net))

network_proteins <- igraph::as_ids(igraph::V(network))
network_proteins <- sort(network_proteins)
return(network_proteins)
}