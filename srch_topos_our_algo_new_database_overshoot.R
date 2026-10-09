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
load("new_initial_info_9_6_phospho_bind.RData")
# saveRDS(c(list(targets),list(topos_to_match)),file="topos_to_match.rds")
tmp<-readRDS("topos_to_match.rds")
topos_to_match <- tmp[[2]]


load("protein_interactions_phospho_bind_7_6_2024.RData")
# 
all_proteins_permut_proteins = all_proteins_permut_proteins[which(!(unlist(all_proteins_permut_values[,1])==unlist(all_proteins_permut_values[,2]))),]
all_proteins_permut_values = all_proteins_permut_values[which(!(unlist(all_proteins_permut_values[,1])==unlist(all_proteins_permut_values[,2]))),]

####################################################################################
############## load separate protein_interaction files and create variables ########

# table_tmp <-   read.csv(file = "protein_interaction_phospho_4_24.csv",
#                         header = TRUE)
# nodes <- sort(unique(c(table_tmp[, 1], table_tmp[, 2])))
# nodes <- unname(sort(nodes))
# protein_range <- 1:length(nodes)
# links<-table_tmp
# 
# net <-
#   igraph::graph_from_data_frame(table_tmp, directed = TRUE, vertices = NULL)
# routes <- list.files("./Strict/Data/only_phospho/", full.names = T)

# load(routes[1])
# paths <- path_proteins
# 
# 
# for (i in 2:length(routes)){
#   
#   load(routes[i])
#   paths <- rbind(paths,path_proteins[2:nrow(path_proteins),])
#   
# }
# proteins_table_final <-
#   as.data.frame(apply(paths[,1:2], 1, function(x)
#     which(unlist(nodes) %in% x)))
# 
# 
# all_proteins_permut_proteins <- cbind(paths[,c(1,2,4)])
# all_proteins_permut_values   <- cbind(paths[,c(1,2,3)])
# 
# 
# all_proteins_permut_values <-
#   cbind(all_proteins_permut_values, rep(0, nrow(all_proteins_permut_values)), rep(0, nrow(all_proteins_permut_values)))
# 
# for (kk in 1:nrow(all_proteins_permut_values)) {
#   
#   
#   if (any(1 %in% unlist(all_proteins_permut_values[kk, 3]))) {
#     all_proteins_permut_values[kk, 4] <- 1
#     
#   }
#   if (any(-1 %in% unlist(all_proteins_permut_values[kk, 3]))) {
#     all_proteins_permut_values[kk, 5] <- -1
#     
#   }
#   
# }
# save(all_proteins_permut_proteins,all_proteins_permut_values,file = "protein_interactions_15_4_2024.RData")

####################################################
# tmp target
# Target_names_exist <- list.files("Results/length_7_phospho_bind_7_6")
# Target_names_exist <- sub("\\.RData$", "", Target_names_exist)
# targets_remained <- targets[!targets %in% Target_names_exist]
# network_targets  <- targets_remained

###############################################################################
# defining variables _________________________________________________________



topos_to_match <- c(0,0,0,-1,0,0,1,1,0)#topos_to_match[, 1:9]
nodes <- unlist(nodes)

# network_targets  <- targets
final_output <- list()

# final_output <- foreach::foreach(k = 1:length(network_targets), .combine = "cbind")%dopar% {

library(foreach)
library(tidyr)



##############################################################################
############### find the topology examples using our algorithm ##########

for (k in 39:length(targets)) { # 5:length(network_targets)
  #
  
  # defining the target
  target <- targets[[k]]
  print(k)
  tryCatch({
    
    if (target%in% nodes) {
      target_node_number <- which(nodes == target)
      output <- list()
      for (i in 1) {
        # browser()
        
        topology      <- topos_to_match
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
             file = paste0("Results/length_7_phospho_bind_7_6_overshoot/",
                           target,
                           ".RData",
                           collapse = ""))
      }
      
      list(output)
      doParallel::stopImplicitCluster()
      gc()
    }
  })
  
}

