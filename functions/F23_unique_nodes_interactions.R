unique_nodes_interactions <-
  function(potential_proteins_table,
           node_range,
           network_proteins,
           num_nodes,
           interaction_matrix,
           all_proteins_permut_proteins,
           all_proteins_permut_values) {
    ############ check if the connection between two nodes is not passed
    ############ through the other nodes
    
    
    bad_rows <- NULL
    
    # a for loop in each row of potential protein tables to see if they are correct
    for (i in 1:nrow(potential_proteins_table)) {
      if (length(unique(as.vector(potential_proteins_table[i,]))) != num_nodes) {
        # browser()
        
        bad_rows <- c(bad_rows, i)
        
      } else{
        # choosing 2 nodes out of X
        node_index_final <- combn(node_range, 2)
        
        
        for (tmp_colname in 1:ncol(node_index_final)) {
          # browser()
          
          w <- node_index_final[1, tmp_colname]
          s <- node_index_final[2, tmp_colname]
          tmp_names <- potential_proteins_table[i, c(node_range[-c(w, s)])]
          
          print(tmp_colname)
          if (interaction_matrix[w, s] != 0) {
            
            
            # browser()
            
            # checking from w to s
            tmp_path_of_proteins_ws <-
              all_proteins_permut_proteins[which(
                all_proteins_permut_proteins[, 1] ==  which(network_proteins == potential_proteins_table[i, w]) &
                  all_proteins_permut_proteins[, 2] ==  which(network_proteins == potential_proteins_table[i, s])
              ), 3][["path_proteins"]]
            
            
            tmp_path_values_ws <-
              all_proteins_permut_values[which(
                all_proteins_permut_proteins[, 1] ==  which(network_proteins == potential_proteins_table[i, w]) &
                  all_proteins_permut_proteins[, 2] ==  which(network_proteins ==
                                                                potential_proteins_table[i, s])
              ), 3][["path_values"]]
            
            
            qq <- 0
            for (tt in 1:length(tmp_path_of_proteins_ws)) {
              # browser()
              
              if ((!any(tmp_names %in% tmp_path_of_proteins_ws[[tt]])) & qq==0) {
                qq <- qq + 1
                
                if (tmp_path_values_ws[[tt]] != interaction_matrix[w, s]) {
                  
                  bad_rows <- c(bad_rows, i)
                }
              }
            }
            if (qq == 0) {
              # browser()
              
              bad_rows <- c(bad_rows, i)
            }
          }
          if (interaction_matrix[s, w] != 0) {
            
            
            # checking from s to w
            tmp_path_of_proteins_sw <-
              all_proteins_permut_proteins[which(
                all_proteins_permut_proteins[, 2] == which(network_proteins == potential_proteins_table[i, w]) &
                  all_proteins_permut_proteins[, 1] ==  which(network_proteins == potential_proteins_table[i, s])
              ), 3][["path_proteins"]]
            
            
            tmp_path_values_sw <-
              all_proteins_permut_values[which(
                all_proteins_permut_proteins[, 2] ==  which(network_proteins == potential_proteins_table[i, w]) &
                  all_proteins_permut_proteins[, 1] ==  which(network_proteins == potential_proteins_table[i, s])
              ), 3][["path_values"]]
            
            qq <- 0
            
            for (tt in 1:length(tmp_path_of_proteins_sw)) {
              if ((!any(tmp_names %in% tmp_path_of_proteins_sw[[tt]])) & qq==0) {
                qq <- qq + 1
                if (tmp_path_values_sw[[tt]] != interaction_matrix[s, w]) {
                  # browser()
                  
                  bad_rows <- c(bad_rows, i)
                }
                
              }
            }
            if (qq == 0) {
              # browser()
              
              bad_rows <- c(bad_rows, i)
            }
          }
        }
      }
    }
    return(bad_rows)
  }
