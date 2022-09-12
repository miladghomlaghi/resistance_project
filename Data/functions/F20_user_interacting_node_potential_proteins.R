user_interacting_nodes_potential_proteins <-
  function(user_node_position,
           node_range,
           interaction_matrix,
           all_proteins_permut_values,
           potential_nodes,
           potential_nodes_links,
           main_table)
    
  {
    ############# extracting the potential networks based on the users inputs #######################
    ###### starting with the user node: checking the potential proteins for the nodes
    # that interact with the user node
    for (w in node_range[-user_node_position]) {
      # w represents other nodes
      
      for (s in user_node_position) {
        # s is user-node number
        
        # finding the nodes that are regulated by the user node
        if (interaction_matrix[s, w] != 0) {
          # find proteins in the table that contains all combination of proteins interaction (all_proteins_permut)
          # to find possible proteins that are regulated by node S in a way that is determined by user (in interaction_matrix)
          # browser()
          
          
          tmp_all_proteins_permut1 <-
            all_proteins_permut_values[which(all_proteins_permut_values[, 1] ==  potential_nodes[[as.character(s)]]), ]
          # browser()
          new_potential_nodes <- list()
          new_potential_nodes <-
            foreach::foreach(kk = 1:nrow(tmp_all_proteins_permut1),
                             .combine = "rbind") %dopar% {
                               tmp <- NULL
                               # for (kk in 1:nrow(tmp_all_proteins_permut1)) {
                                 if (all(interaction_matrix[s, w] %in% unlist(tmp_all_proteins_permut1[kk, 3]))) {
                                   # new_potential_nodes[[length(new_potential_nodes) + 1]]
                                   tmp <-
                                     unlist(tmp_all_proteins_permut1[kk, 2])
                                   
                                 }
                                 tmp
                               
                             }
          
          # check if the chosen interaction by user is not among the network interactions
          # if not then terminate the search
          
          tmp_cond1 <-
            intersect(potential_nodes[[as.character(w)]], unlist(new_potential_nodes))
          
          if (length(tmp_cond1) == 0) {
            # browser()
            return(list("not possible"))
          } else {
            # updating the potential proteins for node W by intersecting the new potential
            # proteins and previous ones
            
            if (w %in% names(main_table)) {
              main_table <-
                main_table[(main_table[[as.character(w)]] %in% tmp_cond1), ]
            }
            
            if (!any(w %in% names(main_table))) {
              column_names <- c(names(main_table), w)
              main_table <-
                expand_grid(main_table, as.data.frame(unique(tmp_cond1)))
              colnames(main_table) <-
                column_names
              
              
            }
          }
        }
        
        # finding the nodes that regulate the user node
        if (interaction_matrix[w, s] != 0) {
          # find proteins in the table that contains all combination of proteins interaction (all_proteins_permut)
          # to find possible proteins that regulate by node S in a way that is determined by user (in interaction_matrix)
          # browser()
          tmp_all_proteins_permut <-
            all_proteins_permut_values[which(all_proteins_permut_values[, 2] ==  potential_nodes[[as.character(s)]]), ]
          
          new_potential_nodes2 <- list()
          
          new_potential_nodes2 <-
            foreach::foreach(kk = 1:nrow(tmp_all_proteins_permut),
                             .combine = "rbind") %dopar% {
                               tmp <- NULL
                               
                               # for (kk in 1:nrow(tmp_all_proteins_permut)) {
                               if (all(interaction_matrix[w, s] %in% unlist(tmp_all_proteins_permut[kk, 3]))) {
                                 # browser()
                                 
                                 # new_potential_nodes2[[length(new_potential_nodes2) + 1]] <-
                                 tmp <-
                                   unlist(tmp_all_proteins_permut[kk, 1])
                                 
                               }
                               tmp
                             }
          
          tmp_cond <-
            intersect(potential_nodes[[as.character(w)]], unlist(new_potential_nodes2))
          
          if (length(tmp_cond) == 0) {
            # browser()
            return(list("not possible"))
          } else {
            if (w %in% names(main_table)) {
              main_table <-
                main_table[(main_table[[as.character(w)]] %in% tmp_cond), ]
            }
            
            if (!any(w %in% names(main_table))) {
              column_names <- c(names(main_table), w)
              main_table <-
                expand_grid(main_table, as.data.frame(unique(tmp_cond)))
              colnames(main_table) <- column_names
              
            }
            
          }
        }
        
      }
      
    }
    return(main_table)
  }