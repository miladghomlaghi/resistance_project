
find_exact_structures_new_4 <-
  function(topology,
           user_node,
           user_node_position,
           nodes,
           links,
           network) {
    #  topology<-c(-1,1,1,0,0,1,0,0,0)
    #   user_node<-"EGFR"
    #    user_node_position <-1
    #   nodes<-nodes[[3]]
    #    links<-links[[3]]
    #   network<-nets[[3]]
    
    source("./data/functions/F02_get_link.R")
    numCores <- bigstatsr::nb_cores()
    doParallel::registerDoParallel(cl <- makeCluster(numCores))
    
    dist <- 5
    positions <- 1:length(topology)
    num_nodes <- sqrt(length(topology))
    node_range <- 1:num_nodes
    network_proteins <- sort(nodes)
    # browser()
    
    # for user_node
    user_protein_number <-
      which(is.element(network_proteins, as.character(user_node)))
    
    node_index_permutation <-
      gtools::permutations(num_nodes, 2, 1:num_nodes, repeats.allowed = T)
    node_index_permutation <- node_index_permutation[, c(2, 1)]
    
    if (length(network_proteins) > 300) {
      near_net <-
        igraph::ego(
          network,
          dist,
          network_proteins = as.character(user_node),
          mode = "out",
          mindist = 0
        )
      network <- igraph::induced_subgraph(network, unlist(near_net))
      
      network_proteins <- igraph::as_ids(igraph::V(network))
      network_proteins <- sort(network_proteins)
    }
    
    ## calculating the existence of a link between any two proteins in this network
    
    all_proteins_permut <-
      gtools::permutations(length(network_proteins),
                           2,
                           1:length(network_proteins),
                           repeats.allowed = F)
    
    # browser()
    
    check_link_exist = 0
    
    for (i in 1:nrow(all_proteins_permut)) {
      check_link_exist[i] <-
        get_link(network, network_proteins[all_proteins_permut[i, 1]], network_proteins[all_proteins_permut[i, 2]], links, full = F)
      
    }
    # browser()
    
    all_proteins_permut = cbind(all_proteins_permut, check_link_exist)
    
    
    ############# extracting the potential networks based on the users inputs #######################
    
    
    #interaction matrix entered by the user
    interaction_matrix = matrix(topology, num_nodes, num_nodes)
    
    #potential proteins for each node
    potential_nodes <- list()
    
    
    # as a start, each node without connection with the users nodes takes all
    #the possible proteins in the network as potential candidates
    for (i in 1:length(user_node_position)){
    potential_nodes[[as.character(user_node_position[[i]])]] <- user_protein_number[[i]]
    }
    for (e in  node_range[-user_node_position]) {
      potential_nodes[[as.character(e)]] <- 1:length(network_proteins)
      
      
    }
    potential_nodes_links <- list()
    
    ###### starting with the user node: checking the potential proteins for the nodes
    # that interact with the user node
    # browser()
    for (w in node_range[-user_node_position]) {
      # w represents other nodes
      
      for (s in user_node_position) {
        # s is user-node number
        
        # finding the nodes that are regulated by the user node
        if (interaction_matrix[s, w] != 0) {
          
          # find proteins in the table that contains all combination of proteins interaction (all_proteins_permut)
          # to find possible proteins that are regulated by node S in a way that is determined by user (in interaction_matrix)
          new_potential_nodes <-
            c(all_proteins_permut[which(all_proteins_permut[, 1] ==  potential_nodes[[as.character(s)]] &
                                all_proteins_permut[, 3] == interaction_matrix[s, w]), 2])
          
          # check if the chosen interaction by user is not among the network interactions
          # if not then terminate the search 
          if (all(new_potential_nodes == 0) |
              (length(potential_nodes[[as.character(w)]]) == 0)) {
            return(list("not possible1"))
          } else {
            # updating the potential proteins for node W by intersecting the new potential
            # proteins and previous ones
            
            potential_nodes[[as.character(w)]] <-
              intersect(potential_nodes[[as.character(w)]], new_potential_nodes)
          }
        }
        
        # finding the nodes that regulate the user node
        if (interaction_matrix[w, s] != 0) {
          
          # find proteins in the table that contains all combination of proteins interaction (all_proteins_permut)
          # to find possible proteins that regulate by node S in a way that is determined by user (in interaction_matrix)
          
          new_potential_nodes2 <-
            c(all_proteins_permut[which(all_proteins_permut[, 2] ==  potential_nodes[[as.character(s)]] &
                                all_proteins_permut[, 3] == interaction_matrix[w, s]), 1])
          
          if (all(new_potential_nodes2 == 0) |
              (length(potential_nodes[[as.character(w)]]) == 0)) {
            # browser()
            return(list("not possible2"))
          } else {
            potential_nodes[[as.character(w)]] <-
              intersect(potential_nodes[[as.character(w)]], new_potential_nodes2)
            
          }
        }
        
        if ((interaction_matrix[w, s] != 0)|
          (interaction_matrix[s, w] != 0)){
        # saving obtained potential links between node S and W 
        nodes_index = paste0(as.character(s), as.character(w))
        
        potential_nodes_links[nodes_index] <-
          list(crossing(var1 = potential_nodes[[as.character(s)]], var2 = potential_nodes[[as.character(w)]]))
        }
      }
      
    }
    
    
    
    # browser()
    ####### In this step we should find potential proteins for nodes that do not have interaction with 
    # the user node. Here we pick these nodes one by one and find potential proteins 
    # based on the proteins that we found in the previous step
    
    # Finding proper nodes to start by looking at the number of potential
    #proteins for each node
    potential_nodes_size <- list()

    for (e in node_range[-user_node_position]) {

      potential_nodes_size[as.character(e)] <-
        length(potential_nodes[[as.character(e)]])
      
    }
    
    loop = order(sapply(potential_nodes_size,"[[",1))
    sequence <- names(potential_nodes_size)[loop]
    
    
    # browser()
    
    #################### starting the search among the non user interactions
   
    for (v in 1:(length(sequence) - 1)) {
      e <- strtoi(sequence[v])
      l <- strtoi(sequence[v+1])
      
      
      a <- NULL
      b <- NULL
      # browser()
      
      if ((interaction_matrix[e, l] != 0) |
          (interaction_matrix[l, e] != 0)) {
        
        combi_nodes1 <-
          crossing(var1 = potential_nodes[[as.character(e)]], var2 = potential_nodes[[as.character(l)]])
        nodes_index = paste0(as.character(e), as.character(l))

        qq <- 0
        
        for (h in 1:nrow(combi_nodes1)) {
          if (combi_nodes1[h, 1] != combi_nodes1[h, 2]) {
            #avoid self loop
            test1 <- FALSE
            test2 <- FALSE
            
            # browser()
            # for link starting from e to l
            test1 <-
              all_proteins_permut[all_proteins_permut[, 1] == combi_nodes1[[h, 1]] &
                          all_proteins_permut[, 2] == combi_nodes1[[h, 2]], 3] == interaction_matrix[e, l]
            
            
            
            # for link starting from l to e
            test2 <-
              all_proteins_permut[all_proteins_permut[, 2] == combi_nodes1[[h, 1]] &
                          all_proteins_permut[, 1] == combi_nodes1[[h, 2]], 3] == interaction_matrix[l, e]
            
            
            
            if (test1 && test2) {
              # browser()
              qq = qq + 1

              potential_nodes_links[[nodes_index]] <-
                combi_nodes1[h, ]
              
              a <- c(a, combi_nodes1[[h, 1]])
              b <- c(b, combi_nodes1[[h, 2]])
              
            }
          }
        }
        if (qq == 0) {
          # browser()
          
          return(list("not possible3"))
        }
        
        # browser()
        
        potential_nodes[[as.character(e)]] <-
          intersect(potential_nodes[[as.character(e)]], a)
        potential_nodes[[as.character(l)]] <-
          intersect(potential_nodes[[as.character(l)]], b)
      }
    }
    
    # browser()
    
    
    ############### now we need to look at the links again and refine them if 
    ############### any protein has been removed from the potential links
  
    
    if (length(loop) > 1) {
      tmp_combinations = combn(loop, 2)
      
      for (v in 1:(ncol(tmp_combinations))) {
        # browser()
        e <- tmp_combinations[1, v]
        l <- tmp_combinations[2, v]
        if ((interaction_matrix[e, l] != 0) |
            (interaction_matrix[l, e] != 0)) {
          nodes_index = paste0(as.character(e), as.character(l))
          
          for (j in (1:nrow(potential_nodes_links[[nodes_index]]))) {
            # browser()
            
            if (!(potential_nodes_links[[nodes_index]][j, 1] %in% potential_nodes[[as.character(e)]]) |
                !(potential_nodes_links[[nodes_index]][j, 2] %in% potential_nodes[[as.character(l)]])) {
              w = potential_nodes_links[[nodes_index]]
              # browser()
              
              potential_nodes_links[[nodes_index]] <- w[-j, ]
            }
          }
        }
      }
      
    }
    
    
    ##############
    ##############
    
    potential_nodes11111 <- list()
    testest <- list()
    for (i in 1:length(potential_nodes_links)) {
      testest[[i]] <- 1:nrow(potential_nodes_links[[i]])
      
    }
    loop_index <- expand.grid(testest)
    
    # browser()
    for (f in 1:nrow(loop_index)) {
      sss <- list()
      
      for (d in 1:length(potential_nodes_links)) {
        # browser()
        
        sss = c(sss, potential_nodes_links[[d]][loop_index[[f, d]],])
      }
      if (length(unique(sss)) == num_nodes) {
        kossher <- NULL
        # browser()
        
        for (d in 1:length(potential_nodes_links)) {
          ppp = substr(names(potential_nodes_links)[d], 1, 1)
          
          if (!(ppp %in% kossher)) {
            potential_nodes11111[[ppp]] <-
              c(potential_nodes11111[[ppp]], potential_nodes_links[[d]][loop_index[[f, d]], 1])
            kossher <- c(kossher, ppp)
          }
          
          ppp = substr(names(potential_nodes_links)[d], 2, 2)
          
          if (!(ppp %in% kossher)) {
            potential_nodes11111[[ppp]] <-
              c(potential_nodes11111[[ppp]], potential_nodes_links[[d]][loop_index[[f, d]], 2])
            
            kossher <- c(kossher, ppp)
          }
        }
        
      }
      
    }
    
    
    ############ check if the connection between two nodes is not passed 
    ############ through the other nodes
    
    
    rr <- list()
    for (i in 1:num_nodes) {
      rr[[i]] = network_proteins[c(unname(unlist(potential_nodes11111[[as.character(i)]])))]
    }

    
    rrr = as.data.frame(rr)
    # browser()
    
    colnames(rrr) <- c(1:ncol(rrr))
    bad <- NULL
    for (i in 1:nrow(rrr)) {
      # browser()
      
      if (length(unique(as.vector(rrr[i, ]))) != num_nodes) {
        bad <- c(bad, i)
        # browser()
        
      } else{
        for (w in node_range) {
          # w represents other nodes
          
          for (s in node_range) {
            # s is user-node number
            # browser()
            
            if (w != s) {
              if (interaction_matrix[w, s] != 0) {
                # finding the nodes that are regulated by the user node

                tmp_names <- rrr[i, c(node_range[-c(w, s)])]
                
                for (h in 1:length(tmp_names)) {
                  tmp_net <- igraph::delete_vertices(network, tmp_names[[h]])
                }
                print(c(i, w, s))
                check_link <-
                  get_link(tmp_net, rrr[i, w], rrr[i, s], links, full =
                             F)
                
                if (check_link != interaction_matrix[w, s]) {
                  bad <- c(bad, i)
                }
              }
            }
          }
        }
      }
    }
    
    # browser()
    
    rrr <- rrr[-bad, ]
    rownames(rrr) <- 1:nrow(rrr)
    
    return(rrr)
    
  }
