




unique_node_initial_search_user_node <-
  function(potential_proteins_table,
           node_range,
           interaction_matrix,
           network_proteins,
           all_proteins_permut_proteins,
           all_proteins_permut_values,
           shorter,
           longer) {

    all_next_suspicious_rows <- NULL
    all_bad_rows <- NULL
    next_suspicious_rows2 <- NULL
    next_suspicious_rows1 <- NULL
    
    #########################################################
    
    if (interaction_matrix[shorter, longer] != 0) {
      # finding index of protein pairs of columns two and three in all_proteins_permut_values variable
      
      index1 <-
        row.match(potential_proteins_table[, c(shorter, longer)], all_proteins_permut_values[, c(1, 2)])
      
      # browser()
      
      tmp_bad_susp <-
        bulk_identify_bad_susp(potential_proteins_table,
                               interaction_matrix,
                               network_proteins,
                               all_proteins_permut_proteins,
                               all_proteins_permut_values,
                               shorter,
                               longer,
                               index1)
      
      
      suspicious_rows <- tmp_bad_susp[[1]]
      all_bad_rows <- tmp_bad_susp[[2]]
      
      
      tmp_susp <- identify_final_susp(
        potential_proteins_table,
        network_proteins,
        all_proteins_permut_proteins[index1[suspicious_rows], 3],
        suspicious_rows,
        index1
      )
      # browser()
      if (size(tmp_susp)[1]>1){
      next_suspicious_rows1 <- unlist(tmp_susp[,1])
      all_bad_rows <- c(all_bad_rows, unlist(tmp_susp[,2]))
      } else{
        next_suspicious_rows1 <- unlist(tmp_susp[1])
        all_bad_rows <- c(all_bad_rows, unlist(tmp_susp[2])) 
      }
      # next_suspicious_rows1 <-  tmp_susp[[1]]
      # all_bad_rows1 <- c(all_bad_rows1, tmp_susp[[2]])
      
      
    }
    
    #######################################################
    
    if (interaction_matrix[longer, shorter] != 0) {
      # finding index of protein pairs in column two and three
      index2 <-
        row.match(potential_proteins_table[, c(longer, shorter)], all_proteins_permut_values[, c(1, 2)])
      
      # browser()
      
      tmp_bad_susp <-
        bulk_identify_bad_susp(
          potential_proteins_table,
          interaction_matrix,
          network_proteins,
          all_proteins_permut_proteins,
          all_proteins_permut_values,
          longer,
          shorter,
          index2
        )
      
      suspicious_rows <- tmp_bad_susp[[1]]
      all_bad_rows <- c(all_bad_rows,tmp_bad_susp[[2]])
      
      
      tmp_susp <- identify_final_susp(
        potential_proteins_table,
        network_proteins,
        all_proteins_permut_proteins[index2[suspicious_rows], 3],
        suspicious_rows
      )
      # browser()
      next_suspicious_rows2 <-  unlist(tmp_susp[,1])
      all_bad_rows <- c(all_bad_rows, unlist(tmp_susp[,2]))
      
      # next_suspicious_rows2 <-  tmp_susp[[1]]
      # all_bad_rows2 <- c(all_bad_rows2, tmp_susp[[2]])
      
    }
      
      ######## remove the redundant rows in all_proteins_permut_proteins to reduce
      ########  the memory usage
      all_next_suspicious_rows <- c(next_suspicious_rows1,next_suspicious_rows2)
      
      
      if (length(all_next_suspicious_rows)!=0){
      index <- index_susp_in_interact_dataset(
        interaction_matrix,
        potential_proteins_table,
        all_proteins_permut_values,
        all_next_suspicious_rows,
        shorter,
        longer
      )
    
      # browser()
      
      
      ############################################################
      ################# check the suspicious rows rigorously
      
      total_bad_rows <-
        unique_node_parallel_search(
          potential_proteins_table[all_next_suspicious_rows,],
          node_range,
          network_proteins,
          interaction_matrix,
          all_proteins_permut_proteins[index, ],
          all_proteins_permut_values[index, ],
          2,
          3
        )
      all_bad_rows <-
        unique(c(all_bad_rows, all_next_suspicious_rows[total_bad_rows]))
      }
    

      
      
      # removing bad rows from the potential proteins table
      
      if (!(is.null(all_bad_rows)|length(all_bad_rows)==0)) {
        potential_proteins_table <-
          potential_proteins_table[-all_bad_rows, ]
      }
      
    
    return(potential_proteins_table)
    
  }


