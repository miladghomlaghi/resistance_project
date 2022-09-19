library(nycflights13)
library(yenpathy)
load("network_data.RData")
load("./Data/topologies/frequency_3nodes.RData")
load("./Data/networks/pathways.RData")
load("./Data/networks/predefined.RData")





# from_node <- sample.int(10,100)
# to_node <- sample.int(10,100)
# 
# small_graph <- data.frame(
#   start = from_node,
#   end = to_node,
#   weight = rep(1,100)
# )
browser()


links <- links[,c("ENTITYA","ENTITYB","EFFECT_BIN")]
links$EFFECT_BIN <-rep(1,nrow(links))
colnames(links)<-c("start","end","weight")

results <-k_shortest_paths(links,
                           from="AKT1",to="IRS1", k = 20)

# ret_links <- rbind(links[["MTOR signaling"]],links[["TGFbeta"]],links[["Hippo"]], links[["EGFR"]],links[["Insulin receptor"]])
ret_links <- rbind(links[["MTOR signaling"]])
ret_links$STATUS <-rep(1,nrow(ret_links))
colnames(ret_links)<-c("start","end","weight")

browser()
ptm <- proc.time()

results <-k_shortest_paths(ret_links,
                 from="IRS1",to="RPS6KB1", k = 1)


####################################################
print(proc.time()-ptm)


ret_links <- rbind(links[["MTOR signaling"]],links[["TGFbeta"]],links[["Hippo"]], links[["EGFR"]],links[["Insulin receptor"]])

# g = graph_from_data_frame(ret_links, directed = F, vertices = unique(c(from_node,to_node)))
graph_pathways<-igraph::graph_from_data_frame(ret_links, directed = FALSE, vertices = NULL)

ptm <- proc.time()

shortest_paths(graph_pathways, from="IRS1",to="TGFBR2", output="vpath")[[1]]

# all_simple_paths(graph_pathways, from="AKT",to="ULK1")

print(proc.time()-ptm)
