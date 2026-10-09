unique_nodes_interactions <-
  function(potential_proteins_table,
           node_range,
           network_proteins,
           interaction_matrix,
           all_proteins_permut_proteins,
           all_proteins_permut_values) {
    ############ check if the connection between two nodes is not passed
    ############ through the other nodes
    
    # steps <- seq(1, nrow(potential_proteins_table), 10000)
    # steps <- c(steps, nrow(potential_proteins_table))
    
    no_cores <- detectCores()
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)
    
    
    
    all_proteins_permut_values <-
      as.data.frame(all_proteins_permut_values)
    
    # browser()
    for (tmp_colname in 2:3) {
      if (nrow(potential_proteins_table) != 0) {
        potential_proteins_table <-
          unique_node_initial_search (
            potential_proteins_table,
            interaction_matrix,
            network_proteins,
            all_proteins_permut_proteins,
            all_proteins_permut_values,
            tmp_colname
          )
      }
    }
    
    # browser()
    
    potential_proteins_table <-
      unique_node_initial_search_user_node (
        potential_proteins_table,
        node_range,
        interaction_matrix,
        network_proteins,
        all_proteins_permut_proteins,
        all_proteins_permut_values,
        2,
        3
      )
    
    
    
    doParallel::stopImplicitCluster()
    gc()
    
    # browser()
    
    
    return(potential_proteins_table)
  }


# save(list=ls(),file = "enviro.rdata")