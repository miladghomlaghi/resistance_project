

unique_node_initial_search <- function(potential_proteins_table,
                                       interaction_matrix,
                                       network_proteins,
                                       all_proteins_permut_proteins,
                                       all_proteins_permut_values,
                                       tmp_colname) {
  node_range <- 1:3
  shorter <- tmp_colname
  longer <- node_range[-c(1, tmp_colname)]
  all_susp_bad_rows<-list(NULL,NULL)
  bad_rows<-NULL
  
  shorter_col_proteins <-
    unlist(unique(potential_proteins_table[, shorter]))
  # browser()

  potential_proteins_table<-potential_proteins_table %>% arrange(.[[shorter]])
  # potential_proteins_table <-
  #   potential_proteins_table[order(potential_proteins_table[, shorter]),]
  

  if (length(shorter_col_proteins) != 0) {
    all_susp_bad_rows <-
      find_bad_susp_initial_search(
        potential_proteins_table,
        interaction_matrix,
        network_proteins,
        all_proteins_permut_proteins,
        all_proteins_permut_values,
        shorter_col_proteins,
        shorter,
        longer
      )
  }
  
  ######## remove the redundant rows in all_proteins_permut_proteins to reduce the memory usage
  
  # bad_rows <- all_susp_bad_rows[[1]]
  # next_suspecious_rows <- all_susp_bad_rows[[2]][[1]]
  # browser()
  if (length(all_susp_bad_rows) > 2) {
    bad_rows <- unique(unlist(all_susp_bad_rows[, 1]))
  } else{
    bad_rows <- unique(unlist(all_susp_bad_rows[1]))
    
  }
  
  
  
  if (length(all_susp_bad_rows) < 3) {
    next_suspecious_rows <- unique(unlist(all_susp_bad_rows[2]))
    
  } else{
    next_suspecious_rows <- unique(unlist(all_susp_bad_rows[, 2]))
  }
  if (length(next_suspecious_rows) != 0) {
    index <- index_susp_in_interact_dataset(
      interaction_matrix,
      potential_proteins_table,
      all_proteins_permut_values,
      next_suspecious_rows,
      shorter,
      1
    )
    
    
    ############################################################
    ################# check the suspicious rows rigorously
    
    total_bad_rows <- NULL
    # browser()
    if (length(index != 0)) {
      total_bad_rows <-
        unique_node_parallel_search(
          potential_proteins_table[next_suspecious_rows,],
          node_range,
          network_proteins,
          interaction_matrix,
          all_proteins_permut_proteins[index, ],
          all_proteins_permut_values[index, ],
          1,
          shorter
        )
      
      # removing bad rows from the potential proteins table
      
      bad_rows <-
        unique(c(next_suspecious_rows[total_bad_rows], bad_rows))
    }
    
  }
  # browser()
  
  if (!(is.null(unlist(bad_rows))))
  {
    potential_proteins_table <-
      potential_proteins_table[-unlist(bad_rows),]
  }
  
  
  # browser()
  
  
  
  return(potential_proteins_table)
}