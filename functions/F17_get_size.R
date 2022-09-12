get_size <- function(network, protA,protB,links_frame,size) {
  
  warning_state <- getOption("warn")
  options(warn = -1)
  
  path<-names(unlist(igraph::shortest_paths(network, protA,protB)[[1]]))
  
  return(length(path)<=size) 
  
  
}
