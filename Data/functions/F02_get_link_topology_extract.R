


get_link_topology_extract <- function(network, nodeA, nodeB, links_frame, full = T) {
  source("./data/functions/F19_calculate_path_value.R")
  
  tmp_path <- list()
  if (nodeA != nodeB) {
    tmp_path <- igraph::all_simple_paths(network, nodeA, nodeB)
    
    path_proteins <- list()
    path_value <- list()
    soring<-list()
    if (length(tmp_path)==0) {
      path_value<-0
      path_proteins<-NULL
      return(c(list(path_value), list(path_proteins)))
      
    } else{
      
      for (i in 1:length(tmp_path)) {
        tmp_path_2 <- names(unlist(tmp_path[[i]]))
        soring[[i]] <- length(tmp_path_2)
        path_proteins[[i]] <- tmp_path_2
        
        path_value[[i]] <-
          calculate_path_value(tmp_path_2, links_frame, full)
        
      }
      # browser()
      
      
      
      length(path_proteins)
      
      
      
      return(c(list(path_value[order(unlist(soring))]), list(path_proteins[order(unlist(soring))])))
    }
  }
}
