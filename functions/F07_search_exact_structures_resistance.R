





find_exact_structures_resistance <-
  function(topology,
           user_node,
           user_node_position,
           nodes,
           links,
           network) {
    # browser()
    
    source("./functions/F02_get_link_topology_extract.R")
    source("./functions/F20_user_interacting_node_potential_proteins.R")
    source("./functions/F21_non_user_interacting_node_potential_proteins.R")
    source("./functions/F22_valid_combinations_test.R")
    source("./functions/F23_unique_nodes_interactions.R")
    
    
    
    numCores <- bigstatsr::nb_cores()
    doParallel::registerDoParallel(cl <- makeCluster(numCores))
    
    dist <- 5
    positions <- 1:length(topology)
    num_nodes <- sqrt(length(topology))
    node_range <- 1:num_nodes
    network_proteins <- sort(nodes)
    
    user_protein_number <- vector()
    
    # for user_node
    for (i in 1:length(user_node)) {
      user_protein_number <- append(user_protein_number,
                                    which(is.element(
                                      network_proteins, as.character(user_node[[i]])
                                    )))
    }
    
    node_index_permutation <-
      gtools::permutations(num_nodes, 2, 1:num_nodes, repeats.allowed = T)
    node_index_permutation <- node_index_permutation[, c(2, 1)]
    
    
    ## calculating the existence of a link between any two proteins in this network
    
    all_proteins_permut <-
      gtools::permutations(length(network_proteins),
                           2,
                           1:length(network_proteins),
                           repeats.allowed = F)
    
    browser()
    check_link_exist = 0
    path_values <- list()
    path_proteins <- list()
    for (i in 1:nrow(all_proteins_permut)) {
      path_values_proteins <-
        get_link_topology_extract(network,
                 network_proteins[all_proteins_permut[i, 1]],
                 network_proteins[all_proteins_permut[i, 2]],
                 links,
                 full = F)
      
      path_values[i] <- list(path_values_proteins[[1]])
      path_proteins[i] <- list(path_values_proteins[[2]])
      
    }
    
    all_proteins_permut_proteins = cbind(all_proteins_permut, path_proteins)
    all_proteins_permut_values = cbind(all_proteins_permut, path_values)
    
    
    ############# extracting the potential networks based on the users inputs #######################
    
    browser()
    
    #interaction matrix entered by the user
    interaction_matrix = matrix(topology, num_nodes, num_nodes)
    
    #potential proteins for each node
    potential_nodes <- list()
    potential_nodes_links <- list()
    
    # as a start, each node without connection with the users nodes takes all
    #the possible proteins in the network as potential candidates
    for (i in 1:length(user_node_position)) {
      # browser()
      potential_nodes[[as.character(user_node_position[[i]])]] <-
        user_protein_number[[i]]
    }
    for (e in  node_range[-user_node_position]) {
      potential_nodes[[as.character(e)]] <- 1:length(network_proteins)
      
    }
    
    if (length(user_node_position) > 1) {
      user_node_comb <- combn(user_node_position, 2)
      # browser()
      
      for (i in 1:ncol(user_node_comb)) {
        if ((interaction_matrix[user_node_comb[[1, i]], user_node_comb[[2, i]]] != 0) |
             (interaction_matrix[user_node_comb[[2, i]], user_node_comb[[1, i]]] != 0)){
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
        potential_nodes_links
      )
    
    if ("not possible1" %in% tmp_user_node_related[[1]]) {
      return(list("not possible1"))
    } else{
      # potential links between user nodes and the other
      potential_nodes_links <- tmp_user_node_related[[1]]
      
      # potential proteins for the user related nodes
      # browser()
      
      potential_nodes <- tmp_user_node_related[[2]]
    }
    
    
    ####### In this step we should find potential proteins for nodes that do not have interaction with
    ####### the user node. Here we pick these nodes one by one and find potential proteins
    ####### based on the proteins that we found in the previous step
    
    
    # Finding proper nodes to start by looking at the number of potential
    #proteins for each node
    
    
    potential_nodes_size <- list()
    
    for (e in node_range[-user_node_position]) {
      potential_nodes_size[as.character(e)] <-
        length(potential_nodes[[as.character(e)]])
      
    }
    
    loop = order(sapply(potential_nodes_size, "[[", 1))
    sequence <- names(potential_nodes_size)[loop]
    # browser()
    
    #################### starting to search among the non user node interactions
    
    tmp_non_user_output <-
      non_user_interacting_nodes_potential_proteins(
        interaction_matrix,
        all_proteins_permut_values,
        potential_nodes,
        potential_nodes_links,
        sequence
      )
    
    if ("not possible3" %in% tmp_non_user_output[[1]]) {
      return(list("not possible3"))
    } else{
      potential_nodes <- tmp_non_user_output$potential_nodes
      potential_nodes_links <-
        tmp_non_user_output$potential_nodes_links
    }
    
    
    # browser()
    
    ############### now we need to look at the links again and refine them if
    ############### any protein has been removed from the potential links
    
    
    potential_nodes_links <- Recheck_links(potential_nodes,
                                           potential_nodes_links,
                                           interaction_matrix,
                                           node_range)
    
    # browser()
    
    ######### 1. choose one of the potential node links for each link 2. see
    ######### the unique number of nodes are equal to the node numbers
    #########if less it means one protein is used in more than one place and if more
    ######### it means one potential node is wrong
    
    
    # Creating a loop in potential node links to see which combination of potential
    # links are valid
    potential_nodes_matching_link <-
      valid_combinations(potential_nodes_links, num_nodes)
    
    ############ check if the connection between two nodes is not passed
    ############ through the other nodes
    
    # create a nice table based on the potential proteins for each node
    tmp_table <- list()
    for (i in 1:num_nodes) {
      tmp_table[[i]] = network_proteins[c(unname(unlist(potential_nodes_matching_link[[as.character(i)]])))]
    }
    potential_proteins_table <- as.data.frame(tmp_table)
    colnames(potential_proteins_table) <-
      c(1:ncol(potential_proteins_table))
    
    # browser()
    bad_rows <-
      unique_nodes_interactions(
        potential_proteins_table,
        node_range,
        network_proteins,
        num_nodes,
        interaction_matrix,
        all_proteins_permut_proteins,
        all_proteins_permut_values
      )
    
    # browser()
    proteins_table_final <- potential_proteins_table
    
    if (!(is.null(bad_rows))) {
      proteins_table_final <-
        proteins_table_final[-unique(bad_rows), ]
      if (nrow(proteins_table_final) != 0) {
        rownames(proteins_table_final) <-
          1:nrow(proteins_table_final)
      }
    }
    return(proteins_table_final)
    
  }
