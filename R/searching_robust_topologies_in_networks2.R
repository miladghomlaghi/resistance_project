




# load("./Strict/Data/results_STAT3")
source("./R/02_load_essential_data.R")
source("./R/01_load_functions.R")
source("./R/00_load_libraries.R")
source("./Strict/F05_get_sign.R")
source("./Strict/F06_get_link.R")
source("./functions/F02_get_link_topology_extract.R")
library(dplyr)
library(foreach)
library(pracma)




# loading the targetable drugs

###############################################################################

# defining variables _________________________________________________________

chalenging_targets <- list()
no_result_targets <- list()
for (k in 1) {
  
  #:length(targets)
  target <- "AKT2"#targets[[i]]
  Max_distance_from_target <- 2
  print(k)
  
  # number of nodes in the topology
  num_nodes <- 3
  # number of edges in the topology
  positions <- 1:num_nodes ^ 2
  node_position <- 1:num_nodes
  
  
  proteins    <- sort(nodes)
  combi_nodes <-
    gtools::permutations(num_nodes, 2, 1:num_nodes, repeats.allowed = T)
  combi_nodes <- combi_nodes[, c(2, 1)]
  range       <- 1:length(proteins)
  
  
  numCores <- bigstatsr::nb_cores()
  doParallel::registerDoParallel(numCores)
  
  # NL_positions <- num_nodes*(p-1)+p
  
  user_node_num <- which(proteins == target)
  # positions <- positions[!is.element(positions,NL_positions)]
  
  # creating a subgraph with maximum distance of 5 from the target node
  sub_net <-
    igraph::ego(
      net,
      Max_distance_from_target,
      nodes = proteins[user_node_num],
      mode = "all",
      mindist = 0
    )
  sub_net <- induced_subgraph(net, unlist(sub_net))
  sub_nodes <- as_ids(V(sub_net))
  
  # browser()
  
  ###############################################################################
  
  ############### find the examples using the original algorithm ################
  
  if (length(sub_nodes) > 2 & length(sub_nodes) < 200) {
    sub_name <- which(sub_nodes == proteins[user_node_num])
    
    # getting all combinations of proteins in the network to use in combination
    # with the target node
    all_combi <-
      gtools::permutations(length(sub_nodes) - 1,
                           num_nodes - 1,
                           range[-sub_name],
                           repeats.allowed = F)
    
    pracma::tic()
    browser()
    # performing parallel computing to test all the possible proteins combination
    # to find the examples matching the topologies
    results <- data.frame()
    god_com <- ""
    
    # results <-
    # foreach::foreach(r = 1:nrow(all_combi), .combine = "rbind") %dopar% {
    for (r in 1:nrow(all_combi)) {
      combi <- rep(0, num_nodes)
      combi[1] <- sub_name
      combi[2] <- all_combi[r, ][1]
      combi[3] <- all_combi[r, ][2]
      combi <- matrix(combi, nrow = 1, ncol = num_nodes)
      combi_names <-
        sub_nodes[combi]
      topo <- rep(0, num_nodes ^ 2)
      # browser()
      for (i in 1:nrow(combi_nodes)) {
        tmp_nodes <- node_position[-combi_nodes[i, ]] #
        tmp_names <- combi_names[tmp_nodes]
        tmp_net <- igraph::delete_vertices(sub_net, tmp_names)
        
        
        
        tmp_link <-
          get_link_topology_extract(tmp_net, combi_names[combi_nodes[i, 1]], combi_names[combi_nodes[i, 2]], links)
        # browser()
        if (!is.null(tmp_link[[1]])) {
          topo[positions[i]] <- unique(tmp_link[[1]][[1]])
        }

      }
      all_topologies <- expand.grid(topo)
      if (ncol(all_topologies)==1){
        if (!is.na( prodlim::row.match(as.list(all_topologies)$Var1, as.data.frame(topos8rob_to_match)))) {
          god_com <- c(god_com, paste0(combi_names, collapse = ","))
        }
        }else{
          browser()
          
      if (!is.na( generics::intersect(as.data.frame(all_topologies), as.data.frame(topos8rob_to_match)))) {
        god_com <- c(god_com, paste0(combi_names, collapse = ","))
      }
        }
      
      data.frame(god_com)
      
    }
  } else{
    chalenging_targets[[i]] <- target
  }
  pracma::toc()
  
  if (exists("results")) {
    results <- as.data.frame(results)
    
    print(i)
    save(
      results,
      file = paste0(
        "./Strict/Data/results_",
        target,
        "_",
        Max_distance_from_target,
        ".RData",
        collapse = ""
      )
    )
  } else {
    # browser()
    print("No results")
    no_result_targets[i] <- target
  }
}
pracma::toc()

###############################################################################

############### find the examples using the new algorithm #####################
# source("./functions/F07_search_exact_structures_resistance.R")
#
# find_exact_structures_resistance(topos8rob_to_match[1,],"STAT3",c(1), sub_nodes, links[,c(1,2,3)], net)

################################################################################
####################### plot the network #######################################
#
#
# #
# subnodes = sub_nodes
# network = sub_net
# links = links
# full = F
#
#
# # browser()
#
#
# subnodes<-unlist(subnodes)
# if(all(subnodes!='')){
#   if(length(subnodes)>150){
#     subnodes = subnodes[sample(1:150,20)]
#     print(subnodes)
#   }
#   sub_net <- get_subnetwork(subnodes,network)
#   sub_net <- igraph::simplify(sub_net, remove.multiple = TRUE, remove.loops = F)
#
#   browser()
#   if(igraph::gsize(sub_net)>0){
#
#
#     net <- visNetwork::toVisNetworkData(sub_net)
#
#     net$edges$group<-links[prodlim::row.match(net$edges,links[,c(1,2)]) ,3]
#     #for (i in 1:nrow(net$edges)){
#     #net$edges$group[i] <- ifelse(net$edges$effect[i])
#     #as.numeric(get_link(sub_net,net$edges[i,1],net$edges[i,2],links))
#
#     # }toVisNetworkData
#
#
#     net$edges$color<-ifelse(net$edges$group==1, "#62adcd","#df7f80")
#     net$edges$arrows.to.type <- ifelse(net$edges$group==1,"arrow", "circle")
#     net$nodes$group<-ifelse(is.element(net$nodes$id,subnodes), "N","B")
#     #
#     if (full){
#       color_all = "#f6f6ee"
#       font_all = '18px arial black'
#     } else {
#       color_all = "#676b77"
#       font_all = '18px arial white'
#     }
#
#
#     visNetwork(nodes = net$nodes, edges = net$edges, height = "1050px",width = "1000px",
#                background = 'white', main = NULL)%>%
#       visEdges(arrows = "to", width = 2, arrowStrikethrough = F)%>%
#       ##7a7d8c"
#       visGroups(groupname = "N", color = list(background = color_all,
#                                               border = "#f6f6ee",  border = "#7a7d8c", highlight = '#e1af28'),
#                 shape="box",shadow = list(enabled = TRUE, size = 25),
#                 font = font_all) %>%
#       visGroups(groupname = "B", color = list(background = "#f6f6ee", border = "#7a7d8c",highlight = '#e1af28'),
#                 shape="box",shadow = list(enabled = TRUE, size = 25),
#                 font = '18px arial black')%>%
#       visNodes(size = 80)%>%
#       visPhysics(solver = "forceAtlas2Based")%>%
#       visInteraction(navigationButtons = TRUE)
#
#   } else{
#
#     sub_net <- visNetwork::toVisNetworkData(sub_net)
#     visNetwork(nodes = sub_net$nodes, edges = sub_net$edges, height = "100px",width = "600px",
#                background = 'white', main = list(text = "No regulations found. Select more proteins.",
#                                                  style = "font-family:Tahoma;
#                                                  color:#822424;font-size:20px;text-align:center;"))
#
#   }}else{
#
#     nodes <- data.frame(id = 1, label="First select at least two molecules in Network's proteins", color='white')
#     edges <- data.frame(from = c(1,1), to = c(1,1), color='white')
#     visNetwork(nodes, edges, height = "100px", width = "600px") %>%
#       visNodes(font = list(size = 22, color='#822424', family="Tahoma"))
#   }