merge_pathways <- function (links, positions){

  
ret_links<-links[[positions[1]]]
  for(i in 1:length(positions)){
    
    ret_links<- rbind(ret_links,links[[positions[i]]])
    
  }
ret_links <-ret_links%>%dplyr::distinct(ENTITYA,ENTITYB,STATUS, .keep_all = TRUE)
nodes_pathways <- unique(c(ret_links$ENTITYA,ret_links$ENTITYB))
graph_pathways <- igraph::graph_from_data_frame(ret_links, directed = TRUE, vertices = NULL)    
net_list <- list(nodes_pathways,ret_links,graph_pathways)
names(net_list)<-c("user_nodes","user_links","user_net")
return(net_list)
}