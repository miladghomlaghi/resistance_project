







find_exact_structures_new <-
  function(topology,
           user_node,
           user_position,
           nodes,
           links,
           network) {
    #  topology<-c(-1,1,1,0,0,1,0,0,0)
    #   user_node<-"EGFR"
    #    user_position <-1
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
    nodes <- sort(nodes)
    
    user_node_num <-
      which(is.element(nodes, as.character(user_node)))
    
    combi_nodes <-
      gtools::permutations(num_nodes, 2, 1:num_nodes, repeats.allowed = T)
    combi_nodes <- combi_nodes[, c(2, 1)]
    
    if (length(nodes) > 300) {
      near_net <-
        igraph::ego(
          network,
          dist,
          nodes = as.character(user_node),
          mode = "out",
          mindist = 0
        )
      network <- igraph::induced_subgraph(network, unlist(near_net))
      
      nodes <- igraph::as_ids(igraph::V(network))
      nodes <- sort(nodes)
    }
    
    ## calculating the existence of a link between any two proteins in this network
    
    all_combi <-
      gtools::permutations(length(nodes) - length(user_node_num),
                           2,
                           1:length(nodes),
                           repeats.allowed = F)
    
    
    
    tmp_link = 0
    
    for (i in 1:nrow(all_combi)) {
      tmp_link[i] <-
        get_link(network, nodes[all_combi[i, 1]], nodes[all_combi[i, 2]], links, full = F)
      
    }
    all_combi = cbind(all_combi, tmp_link)
    
    
    ############# designing the three-node network based on the users input #######################
    
    #interaction matrix entered by the user
    interaction_matrix = matrix(topology, num_nodes, num_nodes)
    
    #potential proteins for each node
    potential_nodes <- list()
    
    
    # as a start, each node without connection with the users nodes takes all
    #the possible proteins in the network as potential candidates
    
    potential_nodes[[as.character(user_position)]] <- user_node_num
    
    for (e in  node_range[-user_position]) {
      potential_nodes[[as.character(e)]] <- 1:length(nodes)
      
      
    }
    potential_nodes_links <- list()
    
    # starting with the user node: checking the possible proteins for the nodes
    # that interact with the user node
    
    for (w in node_range[-user_position]) {
      # w represents other nodes
      
      for (s in user_position) {
        # s is user-node number
        
        if (interaction_matrix[s, w] != 0) {
          # finding the nodes that are regulated by the user node
          
          new_nodes1 <-
            c(all_combi[which(all_combi[, 1] ==  potential_nodes[[as.character(user_position)]] &
                                all_combi[, 3] == interaction_matrix[s, w]), 2])
          
          # if the chosen interaction by user is not among the network interactions
          if (all(new_nodes1 == 0) |
              (length(potential_nodes[[as.character(w)]]) == 0)) {
            return()
          } else{
            potential_nodes[[as.character(w)]] <-
              intersect(potential_nodes[[as.character(w)]], new_nodes1)
          }
        }
        
        
        if (interaction_matrix[w, s] != 0) {
          # finding the nodes that regulate the user node
          
          new_nodes2 <-
            c(all_combi[which(all_combi[, 2] ==  potential_nodes[[as.character(user_position)]] &
                                all_combi[, 3] == interaction_matrix[w, s]), 1])
          
          if (all(new_nodes2 == 0) |
              (length(potential_nodes[[as.character(w)]]) == 0)) {
            return()
          } else {
            potential_nodes[[as.character(w)]] <-
              intersect(potential_nodes[[as.character(w)]], new_nodes2)
            
          }
        }
        index111 = paste0(as.character(s), as.character(w))
        
        potential_nodes_links[index111] <-
          list(crossing(var1 = potential_nodes[[as.character(s)]], var2 = potential_nodes[[as.character(w)]]))
        
      }
      
    }
    
    
    
    
    
    # browser()
    
    # Finding proper nodes to start by looking at the number of potential
    #proteins for each node
    potential_nodes_size <- matrix(length(node_range))
    
    for (e in node_range[-user_position]) {
      # for (l in (e + 1):length(node_range)) {
      
      potential_nodes_size[e] <-
        length(potential_nodes[[as.character(e)]])
      
      
      # }
    }
    
    loop = order(potential_nodes_size)
    
    for (v in 1:(length(loop) - 1)) {
      e <- loop[v]
      l <- loop[v + 1]
      
      
      a <- NULL
      b <- NULL
      if ((interaction_matrix[e, l] != 0) |
          (interaction_matrix[l, e] != 0)) {
        combi_nodes1 <-
          crossing(var1 = potential_nodes[[as.character(e)]], var2 = potential_nodes[[as.character(l)]])
        index111 = paste0(as.character(e), as.character(l))
        # combi_nodes1<-as.data.frame(combi_nodes1)
        
        qq <- 0
        
        for (h in 1:nrow(combi_nodes1)) {
          if (combi_nodes1[h, 1] != combi_nodes1[h, 2]) {
            #avoid self loop
            test1 <- FALSE
            test2 <- FALSE
            
            # browser()
            # for link starting from e to l
            test1 <-
              all_combi[all_combi[, 1] == combi_nodes1[[h, 1]] &
                          all_combi[, 2] == combi_nodes1[[h, 2]], 3] == interaction_matrix[e, l]
            
            
            
            # for link starting from l to e
            test2 <-
              all_combi[all_combi[, 2] == combi_nodes1[[h, 1]] &
                          all_combi[, 1] == combi_nodes1[[h, 2]], 3] == interaction_matrix[l, e]
            
            
            
            if (test1 && test2) {
              # browser()
              qq = qq + 1
              
              potential_nodes_links[[index111]] <-
                combi_nodes1[h,]
              
              a <- c(a, combi_nodes1[[h, 1]])
              b <- c(b, combi_nodes1[[h, 2]])
              
            }
          }
        }
        if (qq == 0) {
          # browser()
          
          return()
        }
        
        
        potential_nodes[[as.character(e)]] <-
          intersect(potential_nodes[[as.character(e)]], a)
        potential_nodes[[as.character(l)]] <-
          intersect(potential_nodes[[as.character(l)]], b)
      }
    }
    
    
    
    # now we need to look at the links again and refine them if any protein has been removed
    if (length(loop) > 1) {
      for (v in 1:(length(loop) - 2)) {
        for (u in 1:(length(loop) - 1)) {
          e <- loop[v]
          l <- loop[u]
          
          index111 = paste0(as.character(e), as.character(l))
          
          for (j in (1:length(potential_nodes_links[[index111]]))) {
            # browser()
            
            if (!(j[1] %in% potential_nodes[[e]]) |
                !(j[2] %in% potential_nodes[[l]])) {
              w = potential_nodes_links[[index111]]
              
              potential_nodes_links[[index111]] <- w[-j]
            }
            
          }
          
        }
      }
      
      
      
      
      potential_nodes11111 <- list()
      for (d in 1:nrow(potential_nodes_links[[1]])) {
        # browser()
        
        for (f in 1:nrow(potential_nodes_links[[2]])) {
          # browser()
          if (length(unique(c(
            potential_nodes_links[[1]][d, ], potential_nodes_links[[2]][f, ]
          ))) == 3) {
            kossher <- NULL
            ppp = substr(names(potential_nodes_links)[1], 1, 1)
            if (!(ppp %in% kossher)) {
              potential_nodes11111[[ppp]] <-
                c(potential_nodes11111[[ppp]], potential_nodes_links[[1]][d, 1])
              kossher <- c(kossher, ppp)
            }
            
            ppp = substr(names(potential_nodes_links)[1], 2, 2)
            if (!(ppp %in% kossher)) {
              potential_nodes11111[[ppp]] <-
                c(potential_nodes11111[[ppp]], potential_nodes_links[[1]][d, 2])
              kossher <- c(kossher, ppp)
            }
            
            ppp = substr(names(potential_nodes_links)[2], 1, 1)
            
            if (!(ppp %in% kossher)) {
              potential_nodes11111[[ppp]] <-
                c(potential_nodes11111[[ppp]], potential_nodes_links[[2]][f, 1])
              kossher <- c(kossher, ppp)
            }
            
            ppp = substr(names(potential_nodes_links)[2], 2, 2)
            
            if (!(ppp %in% kossher)) {
              potential_nodes11111[[ppp]] <-
                c(potential_nodes11111[[ppp]], potential_nodes_links[[2]][f, 2])
              kossher <- c(kossher, ppp)
            }
            
            
          }
          
        }
        
      }
    }
    
    #check if a the connection between two nodes is not passed the other nodes
    rr = list(node1 = nodes[c(unname(unlist(potential_nodes11111$"1")))],
              node2 = nodes[c(unname(unlist(potential_nodes11111$"2")))],
              node3 = nodes[c(unname(unlist(potential_nodes11111$"3")))],
              node4 = nodes[c(unname(unlist(potential_nodes11111$"4")))])
    # browser()
    
    rrr = as.data.frame(rr)
    bad <- NULL
    for (i in 1:nrow(rrr)) {
      for (w in node_range) {
        # w represents other nodes
        
        for (s in node_range) {
          # s is user-node number
          
          if (w != s) {
            if (interaction_matrix[w, s] != 0) {
              # finding the nodes that are regulated by the user node
              
              tmp_names <- rrr[i, node_range[-c(w, s)]]
              
              tmp_net <- igraph::delete_vertices(network, tmp_names)
              print(c(i, w, s))
              tmp_link <-
                get_link(tmp_net, rrr[i, w], rrr[i, s], links, full =
                           F)
              
              # print(c(length(tmp_link),length(interaction_matrix[w, s])))
              if (tmp_link != interaction_matrix[w, s]) {
                bad <- c(bad, i)
              }
              
            }
          }
        }
      }
    }
    # browser()
    
    rrr <- rrr[-bad,]
    rownames(rrr) <- 1:nrow(rrr)
    
    a = 1
    
    return(rrr)
    
  }
