Recheck_links <-
  function(interaction_matrix,
           all_proteins_permut_values,
           user_node_position,
           node_range,
           main_table,
           network_proteins
           
  ) {
    
    
    # Finding proper nodes to start by looking at the number of potential
    # proteins for each node
    
    
    potential_nodes_size <- list()
    # browser()
    
    sequence <- names(lapply(main_table, function(x) length(unique(x))))

    #################### starting to search among the non user node interactions
    
    
    seq_comb <- combn(sequence, 2)
    for (v in 1:ncol(seq_comb)) {
      # browser()
      e <- strtoi(seq_comb[1, v])
      l <- strtoi(seq_comb[2, v])
      
      
      a <- NULL
      b <- NULL
      
      if ((interaction_matrix[e, l] != 0) |
          (interaction_matrix[l, e] != 0)) {
        
        if (e %in% names(main_table)) {
          var1 <- unique(main_table[[as.character(e)]])
          
        } else{
          var1 <- 1:length(network_proteins)
        }
        if (l %in% names(main_table)) {
          var2 <- unique(main_table[[as.character(l)]])
        } else{
          var2 <- 1:length(network_proteins)
        }
        
        combi_nodes1 <-
          expand_grid(var1 = var1, var2 = var2)
        
        qq <- 0
        # browser()
        
        ### if there is no link from one node to the other node we should not investigate
        ### the existence of a link between these two.
        
        
        for (h in 1:nrow(combi_nodes1)) {
          #
          # if (h == 400) {
          #   browser()
          # }
          # browser()
          if (combi_nodes1[h, 1] != combi_nodes1[h, 2]) {
            #avoid self loop
            test1 <- TRUE
            test2 <- TRUE
            
            # for a link starting from e to l check if a the link is valid
            
            if (interaction_matrix[e, l] != 0) {
              link_values1 <-
                all_proteins_permut_values[all_proteins_permut_values[, 1] == combi_nodes1[[h, 1]] &
                                             all_proteins_permut_values[, 2] == combi_nodes1[[h, 2]], 3]
              test1 <-
                (all(interaction_matrix[e, l] %in% link_values1[[1]]))
              
            }
            
            
            # for link starting from l to e
            
            if (interaction_matrix[l, e] != 0) {
              link_values2 <-
                all_proteins_permut_values[all_proteins_permut_values[, 2] == combi_nodes1[[h, 1]] &
                                             all_proteins_permut_values[, 1] == combi_nodes1[[h, 2]], 3]
              test2 <-
                (all(interaction_matrix[l, e] %in% link_values2[[1]]))
              
            }
            nodes_index = paste0(as.character(e), as.character(l))
            
            if (test1 && test2) {
              qq = qq + 1
              # browser()
              
              
              a <- c(a, combi_nodes1[[h, 1]])
              b <- c(b, combi_nodes1[[h, 2]])
              
              
              
            }
          }
        }
        if (qq == 0) {
          return(list("not possible3"))
        }
        # browser()
        
        if(e %in% names(main_table)) {
          main_table <- main_table[(main_table[[as.character(e)]] %in% a),]
        } 
        if(l %in% names(main_table)) {
          main_table <- main_table[(main_table[[as.character(l)]] %in%  b),]
        }

        
      }
    }
    
    return(main_table)
    
    
  }