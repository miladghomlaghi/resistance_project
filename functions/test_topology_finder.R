source('dependencies.R')

functions_names <- list.files("./data/functions", full.names = T)
lapply(functions_names, source)


example_net <- read.csv("network_template.csv", sep = ",", header = F)


load("./data/topologies/topologies_3nodes.RData")
load("./data/topologies/topologies_4nodes.RData")
load("./data/topologies/frequency_3nodes.RData")
load("./data/networks/pathways.RData")
load("./data/networks/predefined.RData")

nodes_file <- unique(c(example_net[, 1], example_net[, 2]))
links_file <- example_net
net_file <-
  igraph::graph_from_data_frame(links_file, vertices = nodes_file, directed = T)

# values$links<-links_file
# values$nodes<-nodes_file
# values$net<-net_file

#######examples for 3-node topology "AMPK"

# topology <- c(0, 0, -1,
#                 0, 0, 0, 
#                 -1, 1, 0)

topology <- c(0, 0, 0,
              0, 0, 0,
              1, 1, 0)
#
user_node <- c("ERK1/2")
user_position <- c(3)


###########examples for 4-node topology "AMPK"

# topology <- c(0, 0, 0, 0,
#               0, 0, 0, 0,
#               1,1, 0, 0,
#               0, 0, -1, 0)
# user_node <- c("AKT")
# user_position <- c(3)
######### for examples for 4-node "MTOR signaling"
# 
# topology <- c(0,-1, 0, 0,
#               1, 0,-1, 0,
#               0, 0, 0,-1,
#               1, 0, 1, 0
#               )
# user_node <- c("AMPK")
# user_position <- c(1)

###### for examples for 5-node "MTOR+EGFR+InsR signaling"
  
# topology <- c(0,-1, 1, 0,-1,
#              -1, 0, 1, 0, 0,
#               0, 0, 0, 0, 0,
#               1, 1, 0, 0, 0,
#               1, 1, 0, 0, 0)
# 
# user_node <- c("IRS1", "EGFR")
# user_position <- c(2, 1)


###### for examples for 5-node "MTOR+EGFR+InsR signaling"

# topology <- c( 0, 0, 0,-1, 0,
#               -1, 0, 0, 0, 0,
#                0,-1, 0, 0, 0,
#                1, 0, 1, 0, 0,
#                0, 0, 0, 1, 0)
# 
# user_node <- c("AKT", "mTORC1")
# user_position <- c(1, 4)

###### for examples for 5-node "MTOR+EGFR+InsR signaling"

# topology <- c( 0, 1,-1, 0, 0,
#                0, 0,-1, 0, 0,
#                0, 1, 0, 0, 0,
#                1, 0, 1, 0, 0,
#                0, 0, 0, 1, 0)
# # 
# user_node <- c("AKT", "SOS1","mTORC1")
# user_position <- c(1, 2, 4)
ptm <- proc.time()

a = find_exact_structures_new_4_alter(topology, user_node, user_position, nodes[["MTOR signaling"]], links[["MTOR signaling"]], nets[["MTOR signaling"]])

# a=find_exact_structures_new_4_alter(topology,user_node,user_position, nodes[["AMPK"]], links[["AMPK"]], nets[["AMPK"]])


## combine two networks
# browser()
# ret_links <- rbind(links[["MTOR signaling"]], links[["EGFR"]],links[["Insulin receptor"]])
# graph_pathways<-igraph::graph_from_data_frame(ret_links, directed = TRUE, vertices = NULL)

# ret_links <- rbind(links[["EGFR"]],links[["Insulin receptor"]])
# graph_pathways<-igraph::graph_from_data_frame(ret_links, directed = TRUE, vertices = NULL)
# # browser()
# # 
# # 
# a = find_exact_structures_new_4_alter(
#   topology,
#   user_node,
#   user_position,
#   # unique(c(c(nodes[["MTOR signaling"]], nodes[["EGFR"]]),nodes[["Insulin receptor"]])),
#   unique(c(nodes[["Insulin receptor"]], nodes[["EGFR"]])),
#     ret_links,
#   graph_pathways)
# print((proc.time() - ptm)/60)

