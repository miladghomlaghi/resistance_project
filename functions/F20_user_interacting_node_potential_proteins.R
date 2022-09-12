user_interacting_nodes_potential_proteins <-
  function(user_node_position,
           node_range,
           interaction_matrix,
           all_proteins_permut_values,
           potential_nodes,potential_nodes_links)
    
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
            all_proteins_permut_values[which(all_proteins_permut_values[, 1] ==  potential_nodes[[as.character(s)]]),]
          
          new_potential_nodes <- list()
          
          for (kk in 1:nrow(tmp_all_proteins_permut1)) {
            if (all(interaction_matrix[s, w] %in% tmp_all_proteins_permut1[kk, 3][[1]])) {
              new_potential_nodes[[length(new_potential_nodes) + 1]] <-
                unlist(tmp_all_proteins_permut1[kk, 2])
              
            }
            
          }
          
          
          # check if the chosen interaction by user is not among the network interactions
          # if not then terminate the search
          
          tmp_cond1 <- intersect(potential_nodes[[as.character(w)]], unlist(new_potential_nodes))
          
          if (length(tmp_cond1)==0) {
            # browser()
            return(list("not possible1"))
          } else {
            # updating the potential proteins for node W by intersecting the new potential
            # proteins and previous ones
            
            potential_nodes[[as.character(w)]] <- tmp_cond1
          }
        }
        
        # finding the nodes that regulate the user node
        if (interaction_matrix[w, s] != 0) {
          # find proteins in the table that contains all combination of proteins interaction (all_proteins_permut)
          # to find possible proteins that regulate by node S in a way that is determined by user (in interaction_matrix)
          # browser()
          tmp_all_proteins_permut <-
            all_proteins_permut_values[which(all_proteins_permut_values[, 2] ==  potential_nodes[[as.character(s)]]),]
          
          new_potential_nodes2 <- list()
          
          for (kk in 1:nrow(tmp_all_proteins_permut)) {
            
            if (all(interaction_matrix[w, s] %in% tmp_all_proteins_permut[kk, 3][[1]])) {
              # browser()
              
              new_potential_nodes2[[length(new_potential_nodes2) + 1]] <-
                unlist(tmp_all_proteins_permut[kk, 1])
              
            }
          }
          
          tmp_cond <- intersect(potential_nodes[[as.character(w)]], unlist(new_potential_nodes2))
          
          if (length(tmp_cond)==0) {
            # browser()
            return(list("not possible1"))
          } else {
            potential_nodes[[as.character(w)]] <- tmp_cond
            
          }
        }
        
        if ((interaction_matrix[w, s] != 0) |
            (interaction_matrix[s, w] != 0)) {
          # saving obtained potential links between node S and W
          nodes_index = paste0(as.character(s), as.character(w))
          
          # browser()
          
          potential_nodes_links[nodes_index] <-
            list(crossing(var1 = potential_nodes[[as.character(s)]], var2 = potential_nodes[[as.character(w)]]))
        }
      }
      
    }
    return(list(potential_nodes_links=potential_nodes_links,potential_nodes=potential_nodes))
  }