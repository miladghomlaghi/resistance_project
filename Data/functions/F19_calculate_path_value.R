calculate_path_value <- function(path,links_frame, full = T) {
  if (is.null(path)) {
    path_value <- 0
  } else{
    path_value <- 1
    for (p in 1:(length(path) - 1)) {
      tmp <-
        links_frame[prodlim::row.match(c(path[p], path[p + 1]), as.data.frame(links_frame[, c(1, 2)])), 3]
      path_value <- c(path_value, tmp)
      
    }
    #path_value<-prod(path_value[-c(1,length(path_value))])
    
    #
    
  }
  
  if (full) {
    return(path_value)
  } else{
    return(prod(path_value))
  }
}
