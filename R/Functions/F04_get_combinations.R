get_combinations <- function (merged_data, topology,nodes_num) {
  

  nodes_combi<- merged_data[merged_data[,1]==topology,3]
  
  if(nodes_num==3){
  nodes_combi_split <- c("Node 1","Node 2", "Node 3")
  }else{
    nodes_combi_split <- c("Node 1","Node 2", "Node 3", "Node 4")
  }
  
  tmp <- nodes_combi
  tmp_split <- unlist(strsplit(tmp,"/"))
  
  for (s in 1: length(tmp_split)){
    tmp2 <- tmp_split[s]
    tmp2_split<-unlist(strsplit(tmp2,","))
    nodes_combi_split <- rbind(nodes_combi_split,tmp2_split)
  }
  
  if(nodes_num==3){
  colnames(nodes_combi_split)<-c("NodeA","NodeB", "NodeC")
  }else{
    colnames(nodes_combi_split)<-c("NodeA","NodeB", "NodeC", "NodeD")
  }
  
  nodes_combi_split<-nodes_combi_split[-1,]
  row.names(nodes_combi_split)<-1:nrow(nodes_combi_split)
  nodes_combi_split <- gsub(";","/",nodes_combi_split)
  nodes_combi_split <- gsub(":",",",nodes_combi_split)
  
  nodes_combi_split <- data.frame(nodes_combi_split)
  nodes_combi_split$NodeA <- as.character(nodes_combi_split$NodeA)
  nodes_combi_split$NodeB <- as.character(nodes_combi_split$NodeB)
  nodes_combi_split$NodeC <- as.character(nodes_combi_split$NodeC)
  if(nodes_num==4){
    nodes_combi_split$NodeD <- as.character(nodes_combi_split$NodeD)
  }
  
  return(nodes_combi_split)
}