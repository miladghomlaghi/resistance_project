Recheck_links <-
  function(potential_nodes,
           potential_nodes_links,
           interaction_matrix,
           node_range) {
    
    ############### now we need to look at the links again and refine them if
    ############### any protein has been removed from the potential links
    
    tmp_node_comb <- combn(node_range,2)

      for (v in 1:(ncol(tmp_node_comb))) {
        
        e <- tmp_node_comb[1, v]
        l <- tmp_node_comb[2, v]
        
        if ((interaction_matrix[e, l] != 0) |
            (interaction_matrix[l, e] != 0)) {
          nodes_index = paste0(as.character(e), as.character(l))
          # browser()
          
          if (nodes_index %in% names(potential_nodes_links))  {
            for (j in (1:nrow(potential_nodes_links[[nodes_index]]))) {
              if (!(potential_nodes_links[[nodes_index]][j, 1] %in% potential_nodes[[as.character(e)]]) |
                  !(potential_nodes_links[[nodes_index]][j, 2] %in% potential_nodes[[as.character(l)]])) {
                w = potential_nodes_links[[nodes_index]]
                # browser()
                
                potential_nodes_links[[nodes_index]] <- w[-j,]
                
              }
            }
          } else{
            nodes_index = paste0(as.character(l), as.character(e))
            for (j in (1:nrow(potential_nodes_links[[nodes_index]]))) {
              if (!(potential_nodes_links[[nodes_index]][j, 1] %in% potential_nodes[[as.character(l)]]) |
                  !(potential_nodes_links[[nodes_index]][j, 2] %in% potential_nodes[[as.character(e)]])) {
                w = potential_nodes_links[[nodes_index]]
                # browser()
                
                potential_nodes_links[[nodes_index]] <- w[-j,]
                
              }
            }
          }

        }
      }
    
    
    
    return(potential_nodes_links)
  }