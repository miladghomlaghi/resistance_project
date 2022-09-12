


# load("./Strict/Data/results_STAT3")
source("./R/02_load_essential_data.R")

source("./R/01_load_functions.R")
source("./R/00_load_libraries.R")
source("./Strict/F05_get_sign.R")
source("./Strict/F06_get_link.R")

library(foreach)
library(pracma)


browser()

# loading the tasrgetable drugs
load(paste0("./Data/matching/", "Targets.RData", collapse = ""))

###############################################################################

# defining variables _________________________________________________________
target <- "STAT3"
Max_distance_from_target <- 5

# number of nodes in the topology
num_nodes <- 3
# number of edges in the topology
positions <- 1:num_nodes^2
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
  igraph::ego(net,
    Max_distance_from_target,
    nodes = proteins[user_node_num],
    mode = "all",
    mindist = 0
  )
sub_net <- induced_subgraph(net, unlist(sub_net))
sub_nodes <- as_ids(V(sub_net))

# browser()


############### find the examples using the original algorithm

if (length(sub_nodes) > 2) {
  sub_name <- which(sub_nodes == proteins[user_node_num])

  # getting all combinations of proteins in the network to use in combination 
  # with the target node
  all_combi <- gtools::permutations(length(sub_nodes) - 1, num_nodes - 1,
    range[-sub_name],
    repeats.allowed = F
  )

  pracma::tic()
  
  # performing parallel computing to test all the possible proteins combination
  # to find the examples matching the topologies

  results <-
    foreach::foreach(r = 1:nrow(all_combi), .combine = "rbind") %dopar% {


      combi <- rep(0, num_nodes)
      combi[1] <- sub_name
      combi[2] <- all_combi[r, ][1]
      combi[3] <- all_combi[r, ][2]
      combi <- matrix(combi, nrow = 1, ncol = num_nodes)
      combi_names <-
        sub_nodes[combi]
      topo <- rep(0, num_nodes^2)

      for (i in 1:nrow(combi_nodes)) {
        tmp_nodes <- node_position[-combi_nodes[i, ]] #
        tmp_names <- combi_names[tmp_nodes]
        tmp_net <- igraph::delete_vertices(sub_net, tmp_names)
        tmp_link <-
          get_link(tmp_net, combi_names[combi_nodes[i, 1]], combi_names[combi_nodes[i, 2]], links)
        if (!is.na(tmp_link)) {
          topo[positions[i]] <- tmp_link
        }
      }

      if (sum(topo != 0) > 2) {
        c(topo, combi_names)
      }
    }
}
pracma::toc()

if (exists("results")) {
  results <- as.data.frame(results)


  save(results,file = paste0("./Strict/Data/results_",target,"_",Max_distance_from_target,".RData",collapse = "")
  )
} else {
  print("No results")
}

#plot the network ########################################################
