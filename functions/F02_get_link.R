get_link <- function(network, nodeA,nodeB,links_frame, full=T) {
  
  # warning_state <- getOption("warn")
  # options(warn = -1)
  # 
  # 
  if(nodeA!=nodeB){
   path<-names(unlist(igraph::shortest_paths(network, nodeA,nodeB)[[1]]))

  if (is.null(path)){
    path_value <- 0
   }else{
     path_value<-1
     for (p in 1:(length(path)-1)){

      tmp <- links_frame[prodlim::row.match(c(path[p],path[p+1]),as.data.frame(links_frame[,c(1,2)])),3]
       path_value <-c(path_value,tmp)
    
     }
     #path_value<-prod(path_value[-c(1,length(path_value))])
     
  # 
   }}else{
     path_value<-links_frame[prodlim::row.match(c(nodeA,nodeB),links_frame[,c(1,2)]),3]
   }
  
  # 
  # options(warn = warning_state)
 # print(path)
 # return(path_value)


  
  if(full){
    return(path_value)
  }else{
    return(prod(path_value))
  }
  
}
