
find_exact_structures <-
  function(topology,
           user_node,
           user_node_position,
           all_proteins_permut_proteins,
           all_proteins_permut_values,
           nodes
           
           ) {
    source("./data/functions/F02_get_link_topology_extract.R")
    source("./data/functions/F20_user_interacting_node_potential_proteins.R")
    source("./data/functions/F21_non_user_interacting_node_potential_proteins.R")
    source("./data/functions/F22_valid_combinations_test.R")
    source("./data/functions/F23_unique_nodes_interactions.R")
    library(parallel)
    
    
    
    # browser()
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
    # browser()
    # if (length(network_proteins) > 300) {
    #   network_proteins <-
    #     big_network_resizing(network_proteins, user_node)
    # }
    
    ## calculating the existence of links between any two proteins in this network
    
    # path_proteins <-
    #   protein_routes(network_proteins, links, network)
    # 
    # all_proteins_permut_proteins <- path_proteins[[1]]
    # all_proteins_permut_values   <- path_proteins[[2]]
    
    
    ############# extracting the potential networks based on the users inputs #######################
    
    
    #interaction matrix entered by the user
    interaction_matrix = matrix(topology, num_nodes, num_nodes)
    # browser()
    
    tmp_user_node_related <- User_interacting_node_protein_finder(
      user_node_position,
      network_proteins,
      user_protein_number,
      node_range,
      interaction_matrix,
      all_proteins_permut_values
    )
    
    if ("not possible" %in% tmp_user_node_related[[1]] | length(tmp_user_node_related[[1]])==0) {
      return(list("not possible"))
    } else{
      # potential links between user nodes and the other
      
      # potential proteins for the user related nodes
      main_table <- tmp_user_node_related
    }
    
    # browser()
    main_table <- main_table[!apply(main_table,1,function(x) any(duplicated(x))),]
    main_table <- main_table[,sort(names(main_table))]
    
    ####### In this step we should find potential proteins for nodes that do not have interaction with
    ####### the user node. Here we pick these nodes one by one and find potential proteins
    ####### based on the proteins that we found in the previous step
    
    
    #################### starting to search among the non user node interactions
    
    tmp_non_user_output <-
      non_user_interacting_nodes_potential_proteins(
        interaction_matrix,
        all_proteins_permut_values,
        user_node_position,
        node_range,
        main_table,
        network_proteins
      )
    
    if ("not possible" %in% tmp_non_user_output[[1]]| length(tmp_user_node_related[[1]])==0) {
      return(list("not possible"))
    } else{
      main_table <- tmp_non_user_output
      
    }
    main_table <- main_table[!apply(main_table,1,function(x) any(duplicated(x))),]
    main_table <- main_table[,sort(names(main_table))]
    # browser()
    # 
    ############### now we need to look at the links again and refine them if
    ############### any protein has been removed from the potential links
    
    print("step 2 done: Non interacting nodes")
    
    # print(main_table)
    
    # main_table <-   Recheck_links(
    #                               interaction_matrix,
    #                               all_proteins_permut_values,
    #                               user_node_position,
    #                               node_range,
    #                               main_table,
    #                               network_proteins
    #                             )
    # browser()
    
    ######### 1. choose one of the potential node links for each link 2. see
    ######### the unique number of nodes are equal to the node numbers
    ######### if less it means one protein is used in more than one place and if more
    ######### it means one potential node is wrong
    print("step 3 done: recheck")
    
    
    
    ############ check if the connection between two nodes is not passed
    ############ through the other nodes


    
    bad_rows <-
      unique_nodes_interactions(
        main_table,
        node_range,
        network_proteins,
        num_nodes,
        interaction_matrix,
        all_proteins_permut_proteins,
        all_proteins_permut_values
      )
    # browser()
    print("step 4 done: Finding bad rows")
    
    proteins_table_final<-as.data.frame(apply(main_table,2,function(x) network_proteins[x]))
    
    if (!(is.null(bad_rows))) {
      proteins_table_final <-
        proteins_table_final[-unique(bad_rows),]
      # browser()
      if (nrow(proteins_table_final) == 0) {
        return("Not possible")
      }else{
        rownames(proteins_table_final) <-
          1:nrow(proteins_table_final)
      }
    }
    # browser()
    
    return(as.data.frame(proteins_table_final))
    
  }
