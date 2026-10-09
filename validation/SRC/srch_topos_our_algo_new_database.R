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
## importing the network
# load("Data/permut.rData")

load("Data/initial_info.RData")
# load("protein_interactions_20_length_27_10_2022.RData")

# all_proteins_permut_proteins = all_proteins_permut_proteins[which(!(all_proteins_permut[,1]==all_proteins_permut[,2])),]
# all_proteins_permut_values = all_proteins_permut_values[which(!(all_proteins_permut[,1]==all_proteins_permut[,2])),]

####################################################################################
############## load separate protein_interaction files and create variables ########

load("permut.rData")
# 
# 
routes <- list.files("./Strict/Data/our_network", full.names = T)

load(routes[1])
paths <- path_proteins


for (i in 2:length(routes)){
  
  load(routes[i])
  paths <- rbind(paths,path_proteins[2:nrow(path_proteins),])
  
}
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
save(all_proteins_permut_proteins,all_proteins_permut_values,file = "protein_interactions_20_length_13_3_2024.RData")
# browser()
# 
# which(nodes=="MTOR")
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
########### Filtering three proteins to test the algorithm
# all_proteins_permut_proteins = 
# all_proteins_permut_proteins[
#   all_proteins_permut_proteins[, 1] == 4 |
#     all_proteins_permut_proteins[, 1] == 777 |
#     all_proteins_permut_proteins[, 1] == 405 |
#     all_proteins_permut_proteins[, 2] == 4 |
#     all_proteins_permut_proteins[, 2] == 777 |
#     all_proteins_permut_proteins[, 2] == 405 
# ,]
# 
# all_proteins_permut_values = 
#   all_proteins_permut_values[
#     all_proteins_permut_proteins[, 1] == 4 |
#       all_proteins_permut_proteins[, 1] == 777 |
#       all_proteins_permut_proteins[, 1] == 405 |
#       all_proteins_permut_proteins[, 2] == 4 |
#       all_proteins_permut_proteins[, 2] == 777 |
#       all_proteins_permut_proteins[, 2] == 405 
#     ,]
# nodes = c("ABL1" ,"PTPN22" ,"HSPB1")
###############################################################################

# defining variables _________________________________________________________



topos_to_match <- topos_to_match[, 1:9]
nodes <- unlist(nodes)

network_targets  <- targets
final_output <- list()

# final_output <- foreach::foreach(k = 1:length(network_targets), .combine = "cbind")%dopar% {

library(foreach)
library(tidyr)



##############################################################################

############### find the topology examples using our algorithm ##########

for (k in 1:length(network_targets)) { # 
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
    output[[i]] <-
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
    save(output,
         file = paste0("Results/length_20_revised/",
                       target,
                       ".RData",
                       collapse = ""))
  }
  
  list(output)
  doParallel::stopImplicitCluster()
  gc()
  })
  
}
