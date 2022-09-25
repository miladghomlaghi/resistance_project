library(igraph)

initial_all_interactions <- read.csv("Data/databases_proteins/phospho_bind_signor_interactions.csv")

initial_all_interactions <- initial_all_interactions[c("ENTITYA","ENTITYB","EFFECT")]

initial_all_interactions$EFFECT[which(startsWith(initial_all_interactions$EFFECT, 'down'))] <- -1
initial_all_interactions$EFFECT[which(startsWith(initial_all_interactions$EFFECT, 'up'))] <- 1
initial_all_interactions$EFFECT <- strtoi(initial_all_interactions$EFFECT)
links <- unique(initial_all_interactions)

nodes <- unique(c(initial_all_interactions$ENTITYA,initial_all_interactions$ENTITYB))

net <- graph_from_data_frame(links, directed = TRUE, vertices = NULL)

save(links,net,nodes, file = "network_info.RData")
