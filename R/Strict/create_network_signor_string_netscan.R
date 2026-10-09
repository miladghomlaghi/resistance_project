#### SIGNOR-STING FULL NETWORK NETSCAN

load("./Strict/Data/SIGNOR_NETSCAN.Rdata")

signor_links<-links
signor_nodes<-nodes
signor_net<-net

load("./Strict/Data/STRING_NETSCAN.Rdata")

string_links<-links
string_nodes<-nodes
string_net<-net

rm(nodes,links, net)


links<- rbind(string_links,signor_links)
links<-links%>% dplyr::distinct(ENTITYA, ENTITYB, .keep_all = TRUE)

nodes<-unique(c(as.character(links$ENTITYA),as.character(links$ENTITYB)))
net <- igraph::graph_from_data_frame(links,vertices = nodes,directed = T)
net<-simplify(net, remove.loops = F)
nodes<- as_ids(V(net))


save(links,nodes,net,file = './Strict/Data/SIRNOR_STRING_NETSCAN.RData')





