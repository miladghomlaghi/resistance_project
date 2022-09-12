





valid_combinations <-
  function(potential_nodes_links, num_nodes) {
    ######### 1. choose one of the potential node links for each link 2. see
    ######### the unique number of nodes are equal to the node numbers
    ######### if less it means one protein is used in more than one place and if more
    ######### it means one potential node is wrong
    
    source("./data/functions/F25_potential_links_permutation.R")
    browser()
    
    
    potential_nodes_matching_link <- list()
    
    # consider all the possible combination among the potential links
    
    tmp <- potential_links_permutation(potential_nodes_links)
    loop_index1<- tmp[[1]]
    potential_nodes_links_unique<- tmp[[2]]
    potential_nodes_links_new<- tmp[[3]]
    len_links<-unlist(lapply(potential_nodes_links_new, nrow))
    
    # browser()
    
    # start to check all the rows in the loop index
    f <- 1
    
    while (f <= (nrow(loop_index1))) {
      
      print(f)
      
      # starting with two links and then increase one by one
      for (selected_links in 2:ncol(loop_index1)) {
        
        # fill with the links that only have one value
        tmp_fth_potential_link_combination <-
          unlist(potential_nodes_links_unique)
        
        
        s <- vector()
        # loop in the potential link elements   
        for (dd in 1:selected_links) {
          tmp_fth_potential_link_combination <-
            c(tmp_fth_potential_link_combination,
              potential_nodes_links_new[[dd]][loop_index1[[f, dd]],])
          
          s <-
            append(s, as.integer(substr(
              names(potential_nodes_links_new)[dd], 1, 1
            )))
          s <-
            append(s, as.integer(substr(
              names(potential_nodes_links_new)[dd], 2, 2
            )))
          
          names_fth_potential_link_combination <-
            c(tmp_fth_potential_link_combination,
              potential_nodes_links_new[[dd]][loop_index1[[f, dd]],])
        }
        
        
        if (length(unique(tmp_fth_potential_link_combination)) > length(unique(s))) {
          # find all rows that start with the nth vector of loop index with higher
          # and remove them
          # browser()
          
          next_value <- loop_index1[f, selected_links] + 1
          if (next_value <= len_links[selected_links]) {
            # browser()
            f <-
              f + which(loop_index1[(f + 1):nrow(loop_index1), selected_links] == next_value)[[1]]
          } else{
            f <- f + 1
          }
          # browser()
          break
        }
      }
      # }
      
      
      # the unique number nodes in the fth combination should be equal to number
      # of nodes to be a valid combination
      
      if (length(unique(tmp_fth_potential_link_combination)) < num_nodes) {
        f <- f + 1
      } else if (length(unique(tmp_fth_potential_link_combination)) == num_nodes) {
        # browser()
        # here we see that this combination is a valid one so we should put the proteins
        # in the related node
        
        ## first we populate tmp_potential_nodes with links with single length
        
        tmp_potential_nodes <- list()
        if (length(potential_nodes_links_unique) != 0) {
          for (selected_links in 1:length(potential_nodes_links_unique)) {
            #pick name of the node from the potential_link variable to use
            # for the potential node indexing
            
            tmp_protein_node_indexing = substr(names(potential_nodes_links_unique)[selected_links], 1, 1)
            
            tmp_potential_nodes[[tmp_protein_node_indexing]] <-
              c(
                tmp_potential_nodes[[tmp_protein_node_indexing]],
                unlist(potential_nodes_links_unique[[selected_links]][1, 1], use.names =
                         FALSE)
              )
            
            tmp_protein_node_indexing = substr(names(potential_nodes_links_unique)[selected_links], 2, 2)
            
            tmp_potential_nodes[[tmp_protein_node_indexing]] <-
              c(
                tmp_potential_nodes[[tmp_protein_node_indexing]],
                unlist(potential_nodes_links_unique[[selected_links]][1, 2], use.names =
                         FALSE)
              )
          }
        }
        
        ## second we populate tmp_potential_nodes with links more than single length
        
        for (selected_links in 1:length(potential_nodes_links_new)) {
          #pick name of the node from the potential_link variable to use
          # for the potential node indexing
          
          tmp_protein_node_indexing = substr(names(potential_nodes_links_new)[selected_links], 1, 1)
          tmp_potential_nodes[[tmp_protein_node_indexing]] <-
            c(tmp_potential_nodes[[tmp_protein_node_indexing]],
              unlist(potential_nodes_links_new[[selected_links]][loop_index1[[f, selected_links]], 1], use.names =
                       FALSE))
          
          tmp_protein_node_indexing = substr(names(potential_nodes_links_new)[selected_links], 2, 2)
          
          tmp_potential_nodes[[tmp_protein_node_indexing]] <-
            c(tmp_potential_nodes[[tmp_protein_node_indexing]],
              unlist(potential_nodes_links_new[[selected_links]][loop_index1[[f, selected_links]], 2], use.names =
                       FALSE))
        }
        
        ## here we check if each node has a unique link
        
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
        f <- f + 1
      }
      
    }
    return(potential_nodes_matching_link)
  }
