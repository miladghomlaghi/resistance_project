

functions_names <- list.files("./data/functions", full.names = T)
lapply(functions_names, source)

# load("./Strict/Data/results_STAT3")
# source("./R/01_load_functions.R")
# source("./R/00_load_libraries.R")
# source("./Strict/F05_get_sign.R")
# source("./Strict/F06_get_link.R")
library(dplyr)
library(foreach)
library(pracma)
library(openxlsx)
library(tidyr)

load("robust_topologies.Rdata")

# browser()

## importing the network

table_tmp <- read.xlsx("Sample_pathway_database_single_file.xlsx")
proteins <- sort(unique(c(table_tmp[, 1], table_tmp[, 2])))
range       <- 1:length(proteins)

links <- table_tmp
network <-
  igraph::graph_from_data_frame(table_tmp, directed = TRUE, vertices = NULL)
# browser()

# loading the targetable drugs
num_nodes <- 3

load(paste0("./Data/matching/", "Targets.RData", collapse = ''))
load(paste0(
  "./Data/topologies/",
  num_nodes,
  "Nodes/topos_all.RData",
  collapse = ''
))
network_targets  <- proteins[is.element(proteins, targets)]



numCores <- bigstatsr::nb_cores()
doParallel::registerDoParallel(numCores)


###############################################################################

# defining variables _________________________________________________________

challenging_targets <- list()
no_result_targets <- list()

## calculating the existence of links between any two proteins in this network
path_proteins <-
  protein_routes(proteins, links, network)

all_proteins_permut_proteins <- path_proteins[[1]]
all_proteins_permut_values   <- path_proteins[[2]]
final_output <- list()

# length(network_targets)
final_output <- foreach::foreach(k = 1:length(network_targets), .combine = "cbind")%dopar% {
  library(foreach)
  library(tidyr)
  
# for (k in 1:length(network_targets) ) { # 

  # defining the target
  target <- network_targets[[k]]
  print(k)

  
  ###############################################################################
  
  ############### find the topology examples using the our algorithm ##########
  
  
  target_node_number <- which(proteins == target)
  output <- list()
  for (i in 1:nrow(robust_topologies[c(1),])) {
    # browser()
    
    topology      <- robust_topologies[i, ]
    user_position <- c(1)
    output[[i]] <-
      find_exact_structures(
        topology,
        target,
        user_position,
        all_proteins_permut_proteins,
        all_proteins_permut_values,
        proteins
      )
    
  }
  # final_output[[network_targets[[k]]]] <- output
  if (length(output)!=0) {
    
    save(
      output,
      file = paste0(
        "Results/",
        target,
        ".RData",
        collapse = ""
      )
    )
  } 
  list(output)
}


################################################################################
################# Saving the results for each targetable node ##################


