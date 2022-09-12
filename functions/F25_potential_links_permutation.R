potential_links_permutation <- function(potential_nodes_links) {
  ### sorting the potential_nodes_links based on the length of the elements
  
  len1 <- unlist(lapply(potential_nodes_links, nrow))
  potential_nodes_links2 <- potential_nodes_links[len1 > 1]
  
  # elements with length one will be hold in another variable because they
  # dont need to be assessed in the loop
  potential_nodes_links_unique <- potential_nodes_links[len1 == 1]
  
  # sorting the order of the potential node links
  len2 <- unlist(lapply(potential_nodes_links2, nrow))
  potential_nodes_links3 <- potential_nodes_links2[order(len2)]
  len3 <- unlist(lapply(potential_nodes_links3, nrow))
  
  # below is just for expand grid
  tmp_potential_nodes_links3 <-
    potential_nodes_links2[order(-len2)]
  
  #finding number of rows of each potential link
  tmp_num_rows_potential_links1 <- list()
  
  for (i in 1:length(tmp_potential_nodes_links3)) {
    tmp_num_rows_potential_links1[[i]] <-
      1:nrow(tmp_potential_nodes_links3[[i]])
  }
  
  loop_index <- expand.grid(tmp_num_rows_potential_links1)
  loop_index1 <- loop_index[, order(ncol(loop_index):1)]
  # browser()
  return(list(loop_index1, potential_nodes_links_unique,potential_nodes_links3))
}