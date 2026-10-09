
plot_subnet2<-function (subnodes, network,links){
  source("./R/Functions/F02_get_sign.R")
  source("./R/Functions/F03_get_link.R")

  sub_net <- get_subnetwork(subnodes,network)
  sub_net <- igraph::simplify(sub_net, remove.multiple = TRUE, remove.loops = FALSE)
  
  # if (!is.vector(subnodes)){
  #   
  #   subnodes<- c(as.character(subnodes[1,1]),
  #                as.character(subnodes[1,2]),
  #                as.character(subnodes[1,3]))
  #   
  # }

 net_viz <- from_adj_matrix(as.matrix(igraph::get.adjacency(sub_net)), mode = "directed")
 render_graph(net_viz)
}