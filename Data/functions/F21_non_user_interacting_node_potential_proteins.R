non_user_interacting_nodes_potential_proteins <-
  function(interaction_matrix,
           all_proteins_permut_values,
           user_node_position,
           node_range,
           main_table,
           network_proteins) {
    # Finding proper nodes to start by looking at the number of potential
    # proteins for each node
    
    
    potential_nodes_size <- list()
    # browser()
    
    sequence <-
      names(lapply(main_table, function(x)
        length(unique(x))))
    for (i in node_range) {
      if (!(i %in% sequence)) {
        sequence <- c(sequence, as.character(i))
      }
    }
    sequence <- sequence[-which(sequence == user_node_position)]
    #################### starting to search among the non user node interactions
    
    
    seq_comb <- combn(sequence, 2)
    for (v in 1:ncol(seq_comb)) {
      e <- strtoi(seq_comb[1, v])
      l <- strtoi(seq_comb[2, v])
      
      
      # a <- NULL
      # b <- NULL
      
      if ((interaction_matrix[e, l] != 0) |
          (interaction_matrix[l, e] != 0)) {
        if (e %in% names(main_table)) {
          var11 <- as.vector(unique(main_table[[as.character(e)]]))
          
        } else{
          var11 <- 1:length(network_proteins)
        }
        if (l %in% names(main_table)) {
          var22 <- as.vector(unique(main_table[[as.character(l)]]))
        } else{
          var22 <- 1:length(network_proteins)
        }
        
        test1 <- 1:nrow(all_proteins_permut_values)
        test2 <- 1:nrow(all_proteins_permut_values)
        
        if (interaction_matrix[e, l] == 1) {
          test1 <-
            which(
              all_proteins_permut_values[, 1] %in%  var11  &
                all_proteins_permut_values[, 2] %in%  var22 &
                all_proteins_permut_values[, 4] == 1
            )
          
        } else if (interaction_matrix[e, l] == -1) {
          # browser()
          
          test1 <-
            which(
              all_proteins_permut_values[, 1] %in%  var11  &
                all_proteins_permut_values[, 2] %in%  var22 &
                all_proteins_permut_values[, 5] == -1
            )
          
        }
        
        if (interaction_matrix[l, e] == 1) {
          test2 <-
            which(
              all_proteins_permut_values[, 2] == var11 &
                all_proteins_permut_values[, 1] == var22 &
                all_proteins_permut_values[, 4] == 1
            )
          
        } else if (interaction_matrix[l, e] == -1) {
          test2 <-
            which(
              all_proteins_permut_values[, 2] == var11 &
                all_proteins_permut_values[, 1] == var22 &
                all_proteins_permut_values[, 5] == -1
            )
          
        }
        
        row_index <- intersect(test1, test2)
        # browser()
        total_comb_result <-
          all_proteins_permut_values[row_index, c(1, 2)]
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        #
        # combi_nodes1 <-
        #   expand.grid(var1 = var11, var2 = var22)
        # 
        # qq <- 0
        # # browser()
        # 
        # ### if there is no link from one node to the other node we should not investigate
        # ### the existence of a link between these two.
        # total_comb_result <-
        #   foreach::foreach(h = 1:nrow(combi_nodes1), .combine = "rbind") %dopar% {
        #     rr <- NULL
        # 
        #     # for (h in 1:nrow(combi_nodes1)) {
        #     #
        #     # if (h == 400) {
        #     #   browser()
        #     # }
        #     # browser()
        #     if (combi_nodes1[h, 1] != combi_nodes1[h, 2]) {
        #       #avoid self loop
        #       test1 <- TRUE
        #       test2 <- TRUE
        # 
        #       # for a link starting from e to l check if a the link is valid
        # 
        #       if (interaction_matrix[e, l] != 0) {
        #         link_values1 <-
        #           all_proteins_permut_values[all_proteins_permut_values[, 1] == combi_nodes1[[h, 1]] &
        #                                        all_proteins_permut_values[, 2] == combi_nodes1[[h, 2]], 3]
        #         test1 <-
        #           (all(interaction_matrix[e, l] %in% link_values1[[1]]))
        # 
        #       }
        # 
        # 
        #       # for link starting from l to e
        # 
        #       if (interaction_matrix[l, e] != 0) {
        #         link_values2 <-
        #           all_proteins_permut_values[all_proteins_permut_values[, 2] == combi_nodes1[[h, 1]] &
        #                                        all_proteins_permut_values[, 1] == combi_nodes1[[h, 2]], 3]
        #         test2 <-
        #           (all(interaction_matrix[l, e] %in% link_values2[[1]]))
        # 
        #       }
        #       nodes_index = paste0(as.character(e), as.character(l))
        # 
        #       if (test1 && test2) {
        #         qq = qq + 1
        #         # browser()
        # 
        # 
        #         # a <- c(a, combi_nodes1[[h, 1]])
        #         # b <- c(b, combi_nodes1[[h, 2]])
        #         rr <- c(combi_nodes1[[h, 1]], combi_nodes1[[h, 2]])
        # 
        # 
        #       }
        #     }
        #     rr
        #   }
        
        if (length(total_comb_result) == 0) {
          return(list("not possible"))
        }
        # browser()
        
        if (e %in% names(main_table)) {
          main_table <-
            main_table[(main_table[[as.character(e)]] %in% total_comb_result[, 1]),]
        }
        if (l %in% names(main_table)) {
          main_table <-
            main_table[(main_table[[as.character(l)]] %in%  total_comb_result[, 2]),]
        }
        if (!any(e %in% names(main_table))) {
          column_names <- c(names(main_table), e)
          main_table <-
            expand_grid(main_table, as.data.frame(unique(total_comb_result[, 1])))
          colnames(main_table) <- column_names
          
        }
        if (!any(l %in% names(main_table))) {
          column_names <- c(names(main_table), l)
          main_table <-
            expand_grid(main_table, as.data.frame(unique(total_comb_result[, 2])))
          colnames(main_table) <- column_names
          
        }
        
        
        
        
        
        
        # potential_nodes[[as.character(e)]] <-
        #   intersect(potential_nodes[[as.character(e)]], a)
        # potential_nodes[[as.character(l)]] <-
        #   intersect(potential_nodes[[as.character(l)]], b)
        
      }
    }
    
    return(main_table)
    
    
  }