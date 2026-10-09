User_interacting_node_protein_finder <-
  function(user_node_position,
           network_proteins,
           user_protein_number,
           node_range,
           interaction_matrix,
           all_proteins_permut_values) {
    #potential proteins for each node
    potential_nodes <- list()
    potential_nodes_links <- list()
    
    
    
    # browser()
    
    
    # as a start, each node without connection with the users nodes takes all
    #the possible proteins in the network as potential candidates
    for (i in 1:length(user_node_position)) {
      potential_nodes[[as.character(user_node_position[[i]])]] <-
        user_protein_number[[i]]
    }
    
    
    for (e in  node_range[-user_node_position]) {
      potential_nodes[[as.character(e)]] <- 1:length(network_proteins)
      
    }
    
    ## creating an initial main table
    main_table <- making_main_interaction_table_initial(interaction_matrix,potential_nodes,node_range,user_node_position)   
    # browser()
    
    
    if (length(user_node_position) > 1) {
      user_node_comb <- combn(user_node_position, 2)
      
      for (i in 1:ncol(user_node_comb)) {
        if ((interaction_matrix[user_node_comb[[1, i]], user_node_comb[[2, i]]] != 0) |
            (interaction_matrix[user_node_comb[[2, i]], user_node_comb[[1, i]]] != 0)) {
          index <-
            paste0(as.character(user_node_comb[[1, i]]),
                   as.character(user_node_comb[[2, i]]))
          
          potential_nodes_links[index] <-
            list(crossing(var1 = user_protein_number[[1]], var2 = user_protein_number[[2]]))
        }
        
      }
    }
    

    tmp_user_node_related <-
      user_interacting_nodes_potential_proteins(
        user_node_position,
        node_range,
        interaction_matrix,
        all_proteins_permut_values,
        potential_nodes,
        potential_nodes_links,
        main_table
      )
    return(tmp_user_node_related)
  }