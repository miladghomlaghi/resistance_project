closeAllConnections()
rm(list = ls())

library(foreach)
library(pracma)
library(igraph)
library(doParallel)
library(openxlsx)

################################################################################
##################### Loading the information

load("Data/Data/minimal_top_thesis.RData")
load("Data/Data/topos_all.RData")
load("Data/Data/frequency_matlab_210220.RData")
# load("./Data/Data/network_thesis.Rdata")
load("Data/initial_info.RData")
load("new_initial_info_9_6_phospho_bind.RData")


topos_to_match <-
  cbind(topos_all[frecb[order(frecb$Mean_BB, decreasing = T), ]$Topology[1:8], ], rep(7, 8),
        frecb[order(frecb$Mean_BB, decreasing = T), ]$Topology[1:8])
names(topos_to_match) <- names(final_good_top)
topos_to_match <- rbind(final_good_top, topos_to_match)
topos_to_match <- topos_to_match[, 1:9]


################################################################################
########## protein of interest


protein_index <- which(nodes == "YAP1")




###############################################################################
#################### creating variables
#
# dist <- 5
# nodes <- sort(unlist(nodes))

pracma::tic()
#
num_nodes <- 3
positions <- 1:num_nodes ^ 2
p <- 1:num_nodes
#
combi_nodes <-
  gtools::permutations(num_nodes, 2, 1:num_nodes, repeats.allowed = T)
combi_nodes <- combi_nodes[, c(2, 1)]


################################################################################
#################### Analysis
user_node_num=which(targets=="CDK4")

for (user_node_num in 45) { #1:length(targets)
  # finding the subgraph
  near_net <-net
  
  near_nodes <- nodes
  range  <- 1:length(near_nodes)
  
  if (length(near_nodes) > 2) {
    near_name <- which(near_nodes == targets[user_node_num])
    
    # all_combi <-
    #   gtools::permutations(length(near_nodes) - 1,
    #                        num_nodes - 1,
    #                        range[-near_name],
    #                        repeats.allowed = F)
    
    all_combi <- expand.grid(range[-near_name],protein_index)
    all_combi<-all_combi[-which(all_combi$Var1==all_combi$Var2),]
    
    numCores <- bigstatsr::nb_cores()
    doParallel::registerDoParallel(numCores)
    # r=1
    results <- foreach::foreach(r=1:nrow(all_combi),.combine = "cbind",.verbose = TRUE) %dopar% {
      
      god_com <- list()
      
      
      # for (r in 1:nrow(all_combi)) {
      print(r)
      source("R/Strict/F05_get_sign.R")
      source("R/Strict/F06_get_link.R")
      
      combi <- rep(0, num_nodes)
      combi[1] <- near_name
      combi[2] <- all_combi[r, ][1]
      combi[3] <- all_combi[r, ][2]
      combi <- matrix(combi, nrow = 1, ncol = num_nodes)
      combi_names <- near_nodes[unlist(combi)]
      topo <- rep(0, num_nodes ^ 2)
      
      for (i in 1:nrow(combi_nodes)) {
        tmp_nodes <- p[-combi_nodes[i, ]]
        tmp_names <- combi_names[tmp_nodes]
        tmp_net <- igraph::delete_vertices(near_net, unlist(tmp_names))
        tmp_link <-
          get_link(tmp_net, unlist(combi_names[combi_nodes[i, 1]]),  unlist(combi_names[combi_nodes[i, 2]]), links)
        if (!is.na(tmp_link)) {
          topo[positions[i]] <- tmp_link
        }
      }
      
      for (ii in 1:nrow(topos_to_match)) {
        
        if (!is.na(prodlim::row.match(topo[as.vector(topos_to_match[ii, ] !=
                                                     0)], as.data.frame(topos_to_match[ii, as.vector(topos_to_match[ii, ] != 0)]), nomatch = NA))) {
          # browser()
          # god_com <- cbind(god_com, as.list(combi_names))
          god_com <- c(as.list(combi_names),ii)
          
        }
      }
      god_com
    }

    unique(results[,4])
    
    
    doParallel::stopImplicitCluster()
    
    if (exists("results")) {
      results <- unique(t(as.data.frame(results)))
      
      colnames(results)<-c("one","two","three","topo")
      
      save(results,file = paste0("results_phospho_only_",near_nodes[near_name],".RData",collapse = ""))
      
    }
  }
}
pracma::toc()

a=NULL
for (ii in 1:nrow(topos_to_match)) {
  
  a[ii]=  prodlim::row.match(topos_to_match[ii,as.vector(topos_to_match[14, ] !=0)],
                     as.data.frame(topos_to_match[14,as.vector( topos_to_match[14, ] !=0)]), nomatch = NA)
}





