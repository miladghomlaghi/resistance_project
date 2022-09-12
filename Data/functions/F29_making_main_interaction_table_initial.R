making_main_interaction_table_initial <-
  function(interaction_matrix,
           potential_nodes,
           node_range,
           user_node_position) {
    # browser()
    all_potential_combinations <- expand.grid(potential_nodes[as.character(user_node_position)])
    
    return(all_potential_combinations)
    
  }