closeAllConnections()
rm(list = ls())
gc()

functions_names <-
  list.files("./Data/functions", full.names = T)
lapply(functions_names, source)

library(dplyr)
library(foreach)
library(pracma)
library(openxlsx)
library(tidyr)
library(ParallelLogger)
library(prodlim)
library(parallel)
library(doParallel)
# importing the network
# load("Data/permut.rData")

load("Data/initial_info.RData")

################################################################################
########## Analysing our network
# 
table_tmp <- read.xlsx("Sample_pathway_database_single_file.xlsx")
nodes <- sort(unique(c(table_tmp[, 1], table_tmp[, 2])))
nodes <- unname(sort(nodes))
protein_range       <- 1:length(nodes)
dist <- 100

links <- table_tmp
net <-
  igraph::graph_from_data_frame(table_tmp, directed = TRUE, vertices = NULL)
# browser()

# loading the targetable drugs
num_nodes <- 3

targets  <- nodes[is.element(nodes, targets)]
targets<- c(targets,"FGFR4")

####################################################################################
############## Calculate the interaction path among the proteins ########

source("extract_k_shortest_function.R")
functions_names <- list.files("./data/functions", full.names = T)
lapply(functions_names, source)


all_proteins_permut <-
  gtools::permutations(length(nodes),
                       2,repeats.allowed = F)

all_proteins_permut<-rbind(all_proteins_permut, cbind(c(1:length(nodes)),
                                                      c(1:length(nodes))))

steps <- seq(1, nrow(all_proteins_permut), 1)

# tmp_all_proteins_permut <-
# (all_proteins_permut[all_proteins_permut[, 1] == 904 , ])

numCores <- bigstatsr::nb_cores()

  doParallel::registerDoParallel(numCores)
  path_proteins <- path_proteins_function(net,
                                          nodes,
                                          links,
                                          all_proteins_permut)
  gc()
  doParallel::stopImplicitCluster()
  
  closeAllConnections()
  









####################################################################################
############## load separate protein_interaction files and create variables ########


paths <- path_proteins


proteins_table_final <-
  as.data.frame(apply(paths[,1:2], 1, function(x)
    which(unlist(nodes) %in% x)))


all_proteins_permut_proteins <- cbind(paths[,c(1,2,4)])
all_proteins_permut_values   <- cbind(paths[,c(1,2,3)])


all_proteins_permut_values <-
  cbind(all_proteins_permut_values, rep(0, nrow(all_proteins_permut_values)), rep(0, nrow(all_proteins_permut_values)))

for (kk in 1:nrow(all_proteins_permut_values)) {


  if (any(1 %in% unlist(all_proteins_permut_values[kk, 3]))) {
    all_proteins_permut_values[kk, 4] <- 1

  }
  if (any(-1 %in% unlist(all_proteins_permut_values[kk, 3]))) {
    all_proteins_permut_values[kk, 5] <- -1

  }

}
# browser()

# which(nodes=="mTORC1")
# node_index <- 1:length(nodes)
# this_protein <- all_proteins_permut_proteins[all_proteins_permut_proteins[,1]==node_index[577],3]
# 
# for (i in 577){
#   print(i)
#   this_protein <- all_proteins_permut_values[all_proteins_permut_values[,1]==node_index[577],3]
# 
#  sum_index[i] <-  sum(!( this_protein %in% 0))
# 
# }

###############################################################################

# defining variables _________________________________________________________



topos_to_match <- topos_to_match[, 1:9]
nodes <- unlist(nodes)

network_targets  <- targets
final_output <- list()



##############################################################################

############### find the topology examples using our algorithm ##########
library(foreach)
library(tidyr)

for (k in  8) { # 1:length(network_targets
  #
  
  # defining the target
  target <- network_targets[[k]]
  print(k)
  tryCatch({
    
    
    target_node_number <- which(nodes == target)
    output <- list()
    for (i in 1) {#:nrow(topos_to_match[36,])
      # browser()
      
      topology      <- topos_to_match[36, ]
      user_position <- c(1)
      output[[k]] <-
        find_exact_structures(
          topology,
          target,
          user_position,
          all_proteins_permut_proteins,
          all_proteins_permut_values,
          nodes
        )
      
    }
    # final_output[[network_targets[[k]]]] <- output
    if (length(output) != 0) {
      # save(output,
      #      file = paste0("Results/length_20_revised/",
      #                    target,
      #                    ".RData",
      #                    collapse = ""))
    }
    
    list(output)
    doParallel::stopImplicitCluster()
    gc()
  })
  
}
file_path <- "output_FGFR4.csv" 
write.csv(output[[8]], file = file_path, row.names = FALSE)
