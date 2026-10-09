closeAllConnections()
rm(list = ls())
gc()


library(foreach)
library(pracma)
library(igraph)
library(doParallel)
library(openxlsx)

################################################################################
##################### Loading the information
table_tmp <-   read.csv(file = "protein_link_phospho_bind_binary_7_6.csv",
                        header = TRUE)
nodes <- sort(unique(c(table_tmp[, 1], table_tmp[, 2])))
nodes <- unname(sort(nodes))
protein_range <- 1:length(nodes)
links <- table_tmp
net <-
  igraph::graph_from_data_frame(table_tmp, directed = TRUE, vertices = NULL)


load(paste0("./Data/matching/", "Targets.RData", collapse = ''))

targets  <- nodes[is.element(nodes, targets)]



###############################################################################
#################### creating variables
#
dist <- 100
nodes <- sort(nodes)


# finding the subgraph

near_nodes <- list()
bad_targets<-list()
for (i in 1:length(targets)) {
  near_net <-
    igraph::ego(net,
                dist,
                nodes = targets[i],
                mode = "out",
                mindist = 0)
  near_net <- igraph::induced_subgraph(net, unlist(near_net))
  
  if (length(igraph::as_ids(V(near_net))) > 2) {
    
    near_nodes <- unique(c(near_nodes, igraph::as_ids(V(near_net))))
    
  } else{
    
    bad_targets <- c(bad_targets,i)
  }
}

targets <- targets[-c(unlist(bad_targets))]
net     <- igraph::induced_subgraph(net, unlist(near_nodes))

links = links[apply(links[, 1:2], 1, function (x)
  all(x %in% near_nodes)),]

nodes <- near_nodes


save(nodes,links,targets,net,file = "new_initial_info_9_6_phospho_bind.RData")
