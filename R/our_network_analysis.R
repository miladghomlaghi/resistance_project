
# load("./Strict/Data/results_STAT3")
source("./R/01_load_functions.R")
source("./R/00_load_libraries.R")
source("./Strict/F05_get_sign.R")
source("./Strict/F06_get_link.R")
source("./functions/F02_get_link_topology_extract.R")
library(dplyr)
library(foreach)
library(pracma)
library(openxlsx)
load("robust_topologies.Rdata")

## importing the network

table_tmp <- read.xlsx("Sample_pathway_database_single_file.xlsx")
proteins <- unique(c(table_tmp[, 1], table_tmp[, 2]))
links <- table_tmp
net <-
  igraph::graph_from_data_frame(table_tmp, directed = TRUE, vertices = NULL)
# browser()

# loading the targetable drugs

data_to_match <- 'SIGNOR_NETSCAN'
num_nodes <- 3

load(paste0("./Data/matching/", "Targets.RData", collapse = ''))
load(paste0(
  "./Data/topologies/",
  num_nodes,
  "Nodes/topos_all.RData",
  collapse = ''
))

network_targets  <- proteins[is.element(proteins, targets)]

###############################################################################

# defining variables _________________________________________________________

challenging_targets <- list()
no_result_targets <- list()

for (k in 1) {
  #:length(targets)
  target <- network_targets[[k]]
  print(k)
  
  # number of edges in the topology
  positions <- 1:num_nodes ^ 2
  node_position <- 1:num_nodes
  
  
  proteins    <- sort(proteins)
  combi_nodes <-
    gtools::permutations(num_nodes, 2, 1:num_nodes, repeats.allowed = T)
  combi_nodes <- combi_nodes[, c(2, 1)]
  range       <- 1:length(proteins)
  
  
  numCores <- bigstatsr::nb_cores()
  doParallel::registerDoParallel(numCores)
  
  browser()
  
  ###############################################################################
  
  ############### find the topology examples using the original algorithm #####
  
  target_node_number <- which(proteins == target)
  
  # getting all combinations of proteins in the network to use in combination
  # with the target node
  all_combi <-
    gtools::permutations(length(proteins) - 1,
                         num_nodes - 1,
                         range[-target_node_number],
                         repeats.allowed = F)
  
  pracma::tic()
  # browser()
  
  # performing parallel computing to test all the possible proteins combination
  # to find the examples matching the topologies
  results <- data.frame()
  god_com <- ""
  
  results <-
  foreach::foreach(r = 1:nrow(all_combi), .combine = "rbind") %dopar% {
  # for (r in 1:nrow(all_combi)) {
    combi <- rep(0, num_nodes)
    combi[1] <- target_node_number
    combi[2] <- all_combi[r,][1]
    combi[3] <- all_combi[r,][2]
    combi <- matrix(combi, nrow = 1, ncol = num_nodes)
    combi_names <-
      proteins[combi]
    topo <- list()
    for (i in 1:nrow(combi_nodes)) {
      
      tmp_nodes <- node_position[-combi_nodes[i,]] #
      tmp_names <- combi_names[tmp_nodes]
      tmp_net <- igraph::delete_vertices(net, tmp_names)
      # browser()
      
      
      
      tmp_link <- 
        get_link_topology_extract(tmp_net, combi_names[combi_nodes[i, 1]], combi_names[combi_nodes[i, 2]], links)
      if (!is.null(tmp_link[[1]])) {
        topo[positions[i]] <- list(unique(unlist(tmp_link[[1]])))
        
      }else{topo[positions[i]]<-0}
      
    }
    all_topologies <- expand.grid(topo)
    
    if (nrow(all_topologies) == 1) {
      for (topo_check in 1:nrow(robust_topologies)) {
        index <- which(robust_topologies[topo_check,] != 0)
        
        if (sum(all_topologies[index] == robust_topologies[topo_check, index]) ==
            length(index)) {
          browser()
          god_com <- c(god_com, paste0(combi_names, collapse = ","))
          
        }
      }
      
    
  } else{
    for (topo_check in 1:nrow(robust_topologies)) {
      index <- which(robust_topologies[topo_check,] != 0)
      
      for (topologies_index in 1:nrow(all_topologies)) {
        
        if (sum(all_topologies[topologies_index,index] == robust_topologies[topo_check, index]) ==length(index)) {
          browser()
          god_com <- c(god_com, paste0(combi_names, collapse = ","))
        }
      }
    }

  }
  
    unique(data.frame(god_com))
}

pracma::toc()

if (exists("results")) {
  results <- as.data.frame(results)
  
  print(k)
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
