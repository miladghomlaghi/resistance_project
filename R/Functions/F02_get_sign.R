get_sign <- function(network, protA,protB,links_frame) {
  path<-names(unlist(igraph::shortest_paths(network, protA,protB)[[1]]))
  if (is.null(path)){
    path_value <- 0
  }else{
    path_value <- 1
    for (p in 1: length(path)){
      tmp <- links_frame[prodlim:row.match(c(path[p],path[p+1]),links_frame[,c(1,2)]),3]
      path_value <-c(path_value,tmp)
    }
  }
  path_value
  return(path_value) 
}