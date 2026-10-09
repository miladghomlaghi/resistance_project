closeAllConnections()
rm(list = ls())

library(foreach)
library(pracma)
library(igraph)
library(doParallel)
library(openxlsx)

################################################################################
##################### Loading the information

load("../../Data/Data/minimal_top_thesis.RData")
load("../../Data/Data/topos_all.RData")
load("../../Data/Data/frequency_matlab_210220.RData")
# load("./Data/Data/network_thesis.Rdata")
load("../../Data/initial_info.RData")
load("../../new_initial_info_9_6_phospho_bind.RData")


topos_to_match <-
  cbind(topos_all[frecb[order(frecb$Mean_BB, decreasing = T), ]$Topology[1:8], ], rep(7, 8),
        frecb[order(frecb$Mean_BB, decreasing = T), ]$Topology[1:8])
names(topos_to_match) <- names(final_good_top)
topos_to_match <- rbind(final_good_top, topos_to_match)
topos_to_match <- topos_to_match[, 1:9]


################################################################################
########## Analysing our network

# table_tmp <-
#   read.csv(file = "../../network_links_phospho_binary.csv",
#            header = TRUE)
# nodes <- sort(unique(c(table_tmp[, 1], table_tmp[, 2])))
# nodes <- sort(nodes)
# 
# dist <- 100
# 
# links <- table_tmp
# net <-
#   igraph::graph_from_data_frame(table_tmp, directed = TRUE, vertices = NULL)
# # browser()
# 
# # loading the targetable drugs
# num_nodes <- 3
# 
targets  <- nodes[is.element(nodes, targets)]
#################################################################################
#################################################################################
# calculating the distance between target and each rebound protein in the
# experimental data

# importing results of rebound analysis of experimental data
rebound_data <-
  read.csv(file = "./bounceback_CDK_shared_two.csv",
           header = TRUE)

protein_distance = data.frame(col1 = character(),
                              col2 = numeric(),
                              stringsAsFactors = FALSE)

protein_name <- NULL
protein_index<- NULL
ind<-NULL

rebound_data_in_network = rebound_data%>%filter(GeneNames %in%nodes )
protein_name_in_netowrk = unique(rebound_data_in_network$GeneNames)
protein_name=rebound_data$GeneNames

for (i in 1:length(protein_name)) {
  if(!(protein_name[i]=="CDK4")){

    protein_index <- c(protein_index,which(nodes == protein_name[i]))
  }else{
    ind<-c(ind,i)
  }
}

protein_name <- unique(protein_name)
protein_index <- unique(protein_index[!is.na(protein_index)])


#################################################################################
#################################################################################
# proteins that do not show rebound and exist in the network
# importing results of rebound analysis of experimental data
# all_data <-
#   read.csv(file = "./pY221121.csv",
#            header = TRUE)
# 
# all_proteins<- unique(all_data[,1])
#   
# no_rebound_proteins<-all_proteins[!(all_proteins%in% rebound_data[,1])]
# 
# protein_name <- NULL
# protein_index<- NULL
# ind<-NULL
# for (i in 1:length(no_rebound_proteins)) {
#   protein_name[i] <- strsplit(no_rebound_proteins[i], '_')[[1]][2]
# }
# for (i in 1:length(protein_name)) {
#   if(!(protein_name[i]=="SRC")){
#     
#     protein_index <- c(protein_index,which(nodes == protein_name[i]))
#   }else{
#     ind<-c(ind,i)
#   }
# }
# protein_name <- unique(protein_name)
# protein_index <- unique(protein_index[!is.na(protein_index)])

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
target=which(targets=="CDK4")
which(nodes=="YAP1")
# near name: target index in network
# near_nodes: all nodes in the network


for (user_node_num in target) { #1:length(targets)
  # finding the subgraph
  near_net <-net
  
  near_nodes <- nodes
  range  <- 1:length(near_nodes)
  
  if (length(near_nodes) > 2) {
    near_name <- which(near_nodes == targets[user_node_num][[1]])
    
    # all_combi <-
    #   gtools::permutations(length(near_nodes) - 1,
    #                        num_nodes - 1,
    #                        range[-near_name],
    #                        repeats.allowed = F)
    # all_combi <- expand.grid(range[-near_name],protein_index[72])
    all_combi <- expand.grid(range[-near_name],protein_index)
    all_combi<-all_combi[-which(all_combi$Var1==all_combi$Var2),]
    
    numCores <- bigstatsr::nb_cores()
    doParallel::registerDoParallel(numCores)
    # r=1
    results <- foreach::foreach(r=1:nrow(all_combi),.combine = "cbind",.verbose = TRUE) %dopar% {
      
      god_com <- list()
      
      
      # for (r in 1:nrow(all_combi)) {
      print(r)
      source("../../R/Strict/F05_get_sign.R")
      source("../../R/Strict/F06_get_link.R")
      
      combi <- rep(0, num_nodes)
      combi[1] <- near_name
      combi[2] <- all_combi[r, ][1]
      combi[3] <- all_combi[r, ][2]
      combi <- matrix(combi, nrow = 1, ncol = num_nodes)
      combi_names <- near_nodes[unlist(combi)]
      topo <- rep(0, num_nodes ^ 2)
      
      for (i in 1:nrow(combi_nodes)) {
        tmp_nodes <- p[-combi_nodes[i, ]]
        tmp_names <- unlist(combi_names[tmp_nodes])
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
    
    
    
    doParallel::stopImplicitCluster()
    
    if (exists("results")) {
      results <- unique(t(as.data.frame(results)))
      
      colnames(results)<-c("one","two","three","topo")
      
      save(results,file = paste0("results_phospho_only_",near_nodes[near_name],".RData",collapse = ""))
      
    }
  }
}
pracma::toc()



results<-as.data.frame(results)
saveRDS(c(results,protein_name,protein_index),"CDK_result_non_rebound.rds")
 

high_rob_proteins<- (unique(results[(results$topo>=35),c("three","topo")]))
Unique_node_3 <- unique(results[,c("three","topo")])
not_high_rob_proteins<-(Unique_node_3[!(Unique_node_3$three%in% high_rob_proteins$three),])

table(unlist(not_high_rob_proteins$three))

unique(Unique_node_3$three)

unique(results$three)%in%protein_name

unique(high_rob_proteins$three)

nrow(distinct(high_rob_proteins))
length(unique(not_high_rob_proteins$three))
length(unique(results$three))
protein_name_in_netowrk[!protein_name_in_netowrk %in% unique(results$three)]

