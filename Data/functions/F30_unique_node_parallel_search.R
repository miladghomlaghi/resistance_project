unique_node_parallel_search <- function(potential_proteins_table,
                                        node_range,
                                        network_proteins,
                                        interaction_matrix,
                                        all_proteins_permut_proteins,
                                        all_proteins_permut_values,col1,col2) {
  
  
  ############## defining parallel computation

  
  # a for loop in each row of potential protein tables to see if they are correct
  # browser()
  total_bad_rows <-
  foreach::foreach(i = 1:nrow(potential_proteins_table), .combine = "rbind") %dopar% {
  # for (i in 1:nrow(potential_proteins_table)) { #
      
      bad_rows <- NULL

      
        # browser()
        
        w <- col1
        s <- col2
        
        tmp_names <-
          network_proteins[unlist(potential_proteins_table[i, c(node_range[-c(w, s)])])]
        
        W_protein <- potential_proteins_table[[i, w]]
        S_protein <- potential_proteins_table[[i, s]]
        
        
        if (interaction_matrix[w, s] != 0) {
          
          # checking from w to s
          tmp_path_of_proteins_ws <-
            all_proteins_permut_proteins[which(
              all_proteins_permut_proteins[, 1] == W_protein &
                all_proteins_permut_proteins[, 2] ==  S_protein
            ), 3]
          
          
          tmp_path_values_ws <-
            all_proteins_permut_values[which(
              all_proteins_permut_proteins[, 1] ==  W_protein &
                all_proteins_permut_proteins[, 2] ==  S_protein
            ), 3]
          
          
          qq <- 0
          # browser()
          
          for (tt in 1:length(tmp_path_of_proteins_ws[[1]])) {
            # browser()
            
            if ((!any(tmp_names %in% tmp_path_of_proteins_ws[[1]][[tt]])) &
                qq == 0) {
              qq <- qq + 1
              
              if (tmp_path_values_ws[[1]][[tt]] != interaction_matrix[w, s]) {
                # browser()
                
                bad_rows <-  i
              }
              break
            }
          }
          if (qq == 0) {
            # browser()
            
            bad_rows <-  i
          }
        }
        if (interaction_matrix[s, w] != 0) {
          # browser()
          
          # checking from s to w
          tmp_path_of_proteins_sw <-
            all_proteins_permut_proteins[which(
              all_proteins_permut_proteins[, 2] == W_protein &
                all_proteins_permut_proteins[, 1] == S_protein
            ), 3]
          
          
          tmp_path_values_sw <-
            all_proteins_permut_values[which(
              all_proteins_permut_proteins[, 2] ==  W_protein &
                all_proteins_permut_proteins[, 1] ==  S_protein
            ), 3]
          
          qq <- 0
          # browser()
          for (tt in 1:length(tmp_path_of_proteins_sw[[1]])) {
            if ((!any(tmp_names %in% tmp_path_of_proteins_sw[[1]][[tt]])) &
                qq == 0) {
              qq <- qq + 1
              
              if (any(tmp_path_values_sw[[1]][[tt]] != interaction_matrix[s, w])) {
                # browser()
                
                bad_rows <- i
              }
              break
            }
          }
          if (qq == 0) {
            # browser()
            
            bad_rows <- i
          }
        }
      # }
      
      bad_rows
    }

  # closeAllConnections()
  
  return(total_bad_rows)
  
}