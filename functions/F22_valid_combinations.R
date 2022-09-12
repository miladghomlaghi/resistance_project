




valid_combinations <-
  function(potential_nodes_links, num_nodes) {
    ######### 1. choose one of the potential node links for each link 2. see
    ######### the unique number of nodes are equal to the node numbers
    ######### if less it means one protein is used in more than one place and if more
    ######### it means one potential node is wrong
    
    
    potential_nodes_matching_link <- list()
    #finding number of rows of each potential link
    tmp_num_rows_potential_links <- list()
    for (i in 1:length(potential_nodes_links)) {
      tmp_num_rows_potential_links[[i]] <-
        1:nrow(potential_nodes_links[[i]])
      
    }
    loop_index <- expand.grid(tmp_num_rows_potential_links)
    
    # browser()
    for (f in 1:nrow(loop_index)) {
      #f is fth possible combination of potential links
      
      tmp_fth_potential_link_combination <- list()
      
      for (d in 1:length(potential_nodes_links)) {
        #d is dth potential node link
        
        # picking the fth possible combination of potential links for investigation
        tmp_fth_potential_link_combination = c(tmp_fth_potential_link_combination,
                                               potential_nodes_links[[d]][loop_index[[f, d]],])
      }
      
      # the unique number nodes in the fth combination should be equal to number
      # of nodes to be a valid combination
      
      if (length(unique(tmp_fth_potential_link_combination)) == num_nodes) {
        # browser()
        # here we see that this combination is a valid one so we should put the proteins
        # in the related node
        
        # check if a protein is not assigned twice to a node
        assign_check <- NULL
        tmp_potential_nodes <- list()
        for (d in 1:length(potential_nodes_links)) {
          #pick name of the node from the potential_link variable to use
          # for the potential node indexing
          
          tmp_protein_node_indexing = substr(names(potential_nodes_links)[d], 1, 1)
          
          tmp_potential_nodes[[tmp_protein_node_indexing]] <-
            c(tmp_potential_nodes[[tmp_protein_node_indexing]],
              unlist(potential_nodes_links[[d]][loop_index[[f, d]], 1], use.names =
                       FALSE))
          
          tmp_protein_node_indexing = substr(names(potential_nodes_links)[d], 2, 2)
          
          tmp_potential_nodes[[tmp_protein_node_indexing]] <-
            c(tmp_potential_nodes[[tmp_protein_node_indexing]],
              unlist(potential_nodes_links[[d]][loop_index[[f, d]], 2], use.names =
                       FALSE))
        }
        
        tmp_unique_nodes <- lapply(tmp_potential_nodes, unique)
        
        if (!any(lapply(tmp_unique_nodes, length) > 1)) {
          # browser()
          
          if (length(potential_nodes_matching_link) == 0) {
            potential_nodes_matching_link <- tmp_unique_nodes
            
          } else{
            # browser()
            
            potential_nodes_matching_link <-
              Map(c,
                  potential_nodes_matching_link,
                  tmp_unique_nodes)
          }
        }
        
        
        
        # if (!(tmp_protein_node_indexing %in% assign_check)) {
        #   # update potential_nodes_matching_link variable by adding the successful protein to the node variable
        #
        #   potential_nodes_matching_link[[tmp_protein_node_indexing]] <-
        #     c(potential_nodes_matching_link[[tmp_protein_node_indexing]],
        #       potential_nodes_links[[d]][loop_index[[f, d]], 1])
        #
        #   # update assign_check variable to show the related node has been filled in this round
        #   assign_check <-
        #     c(assign_check, tmp_protein_node_indexing)
        # }
        #
        # tmp_protein_node_indexing = substr(names(potential_nodes_links)[d], 2, 2)
        #
        # if (!(tmp_protein_node_indexing %in% assign_check)) {
        #   potential_nodes_matching_link[[tmp_protein_node_indexing]] <-
        #     c(potential_nodes_matching_link[[tmp_protein_node_indexing]],
        #       potential_nodes_links[[d]][loop_index[[f, d]], 2])
        #
        #   assign_check <-
        #     c(assign_check, tmp_protein_node_indexing)
        # }
      }
      
    }
    
    
    return(potential_nodes_matching_link)
  }
