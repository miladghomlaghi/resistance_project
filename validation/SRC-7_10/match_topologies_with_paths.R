closeAllConnections()
rm(list = ls())

library(foreach)
library(pracma)
library(igraph)
library(doParallel)
library(openxlsx)
library(dplyr)

################################################################################
##################### Loading the information

load("../../Data/Data/minimal_top_thesis.RData")
load("../../Data/Data/topos_all.RData")
load("../../Data/Data/frequency_matlab_210220.RData")
# load("./Data/Data/network_thesis.Rdata")
load("../../Data/initial_info.RData")

load("new_initial_info_9_6_phospho_bind.RData")

topos_to_match <-
  cbind(topos_all[frecb[order(frecb$Mean_BB, decreasing = T), ]$Topology[1:8], ], rep(7, 8),
        frecb[order(frecb$Mean_BB, decreasing = T), ]$Topology[1:8])
names(topos_to_match) <- names(final_good_top)
topos_to_match <- rbind(final_good_top, topos_to_match)
topos_to_match <- topos_to_match[, 1:9]

## topo_id >= 35 are the 8 highly robust topologies; 1..34 are the minimal set.

################################################################################
########## Analysing our network

targets  <- nodes[is.element(nodes, targets)]

#################################################################################
# importing results of rebound analysis of experimental data

rebound_data <-
  read.csv(file = "./bounceback_src_shared_two.csv",
           header = TRUE)

protein_name <- NULL
protein_index <- NULL
ind <- NULL
for (i in 1:dim(rebound_data)[1]) {
  protein_name[i] <- strsplit(rebound_data[i, 1], '_')[[1]][1]
}
for (i in 1:length(protein_name)) {
  if (!(protein_name[i] == "SRC")) {
    protein_index <- c(protein_index, which(nodes == protein_name[i]))
  } else {
    ind <- c(ind, i)
  }
}
protein_name  <- unique(protein_name)
protein_index <- unique(protein_index[!is.na(protein_index)])

###############################################################################
#################### creating variables

pracma::tic()

num_nodes <- 3
positions <- 1:num_nodes ^ 2
p <- 1:num_nodes

combi_nodes <-
  gtools::permutations(num_nodes, 2, 1:num_nodes, repeats.allowed = T)
combi_nodes <- combi_nodes[, c(2, 1)]

################################################################################
#################### one-iteration worker as a NAMED FUNCTION
## Returns either a data.frame (>=1 row) or NULL. No rbind anywhere inside.
## Defined at top level so it can be reused by the serial debug harness below.

match_one <- function(r, all_combi, near_name, near_nodes, near_net,
                      combi_nodes, positions, topos_to_match,
                      num_nodes, p, links) {
  
  combi <- rep(0, num_nodes)
  combi[1] <- near_name
  combi[2] <- all_combi[r, ][[1]]
  combi[3] <- all_combi[r, ][[2]]
  combi <- matrix(combi, nrow = 1, ncol = num_nodes)
  ## as.character() returns a PLAIN unnamed character vector (proven in debug).
  cn <- as.character(near_nodes[unlist(combi)])
  
  topo  <- rep(0, num_nodes ^ 2)
  paths <- vector("list", num_nodes ^ 2)
  
  for (i in 1:nrow(combi_nodes)) {
    tmp_nodes <- p[-combi_nodes[i, ]]
    tmp_names <- cn[tmp_nodes]
    tmp_net   <- igraph::delete_vertices(near_net, tmp_names)
    tmp_res   <- get_link_path(tmp_net, cn[combi_nodes[i, 1]],
                               cn[combi_nodes[i, 2]], links)
    if (!is.na(tmp_res$sign)) {
      topo[positions[i]]    <- tmp_res$sign
      paths[[positions[i]]] <- tmp_res$path
    }
  }
  
  row_list <- list()
  for (ii in 1:nrow(topos_to_match)) {
    nz <- as.vector(topos_to_match[ii, ] != 0)
    if (!is.na(prodlim::row.match(topo[nz],
                                  as.data.frame(topos_to_match[ii, nz]), nomatch = NA))) {
      
      for (pos in which(nz)) {
        fr  <- combi_nodes[pos, 1]
        to  <- combi_nodes[pos, 2]
        pth <- paths[[pos]]
        row_list[[length(row_list) + 1]] <- data.frame(
          X1 = cn[1], X2 = cn[2], X3 = cn[3],
          topo_id = as.integer(ii),
          highly_robust = ii >= 35,
          edge_pos = as.integer(pos),
          from_slot = as.integer(fr), to_slot = as.integer(to),
          from_name = cn[fr], to_name = cn[to],
          sign = as.numeric(topo[pos]),
          path = if (length(pth) >= 1) paste(pth, collapse = ">") else "",
          path_len = if (length(pth) >= 1) as.integer(length(pth) - 1) else NA_integer_,
          n_intermediate = as.integer(max(0, length(pth) - 2)),
          intermediates = if (length(pth) > 2)
            paste(pth[2:(length(pth) - 1)], collapse = ">") else "",
          stringsAsFactors = FALSE)
      }
    }
  }
  
  if (length(row_list) == 0) return(NULL)
  dplyr::bind_rows(row_list)   # not rbind -> cannot throw match.names
}

################################################################################
#################### Analysis

target <- which(targets == "SRC")

for (user_node_num in target) {
  
  near_net   <- net
  near_nodes <- nodes
  range      <- 1:length(near_nodes)
  
  if (length(near_nodes) > 2) {
    near_name <- which(near_nodes == targets[user_node_num][[1]])
    
    all_combi <- expand.grid(range[-near_name], protein_index)
    all_combi <- all_combi[-which(all_combi$Var1 == all_combi$Var2), ]
    
    numCores <- bigstatsr::nb_cores()
    doParallel::registerDoParallel(numCores)
    
    ## collect each worker's result into a LIST (combine = c), bind ONCE after.
    res_list <- foreach::foreach(
      r = 1:nrow(all_combi),
      .combine = "c",
      .multicombine = TRUE,
      .packages = c("igraph", "prodlim", "dplyr"),
      .verbose = TRUE) %dopar% {
        
        source("../../R/Strict/F05_get_sign.R")
        source("../../R/Strict/F06_get_link.R")
        source("../../R/Strict/F07_get_link_path.R")
        
        out <- match_one(r, all_combi, near_name, near_nodes, near_net,
                         combi_nodes, positions, topos_to_match,
                         num_nodes, p, links)
        list(out)            # wrap so combine="c" builds a list of frames/NULLs
      }
    
    doParallel::stopImplicitCluster()
    
    ## drop NULLs, then a single safe bind
    res_list <- res_list[!vapply(res_list, is.null, logical(1))]
    results  <- if (length(res_list)) dplyr::bind_rows(res_list) else NULL
    
    if (!is.null(results) && nrow(results) > 0) {
      results <- as.data.frame(results) %>% dplyr::distinct()
      save(results,
           file = paste0("results_phospho_paths_", near_nodes[near_name], ".RData",
                         collapse = ""))
    }
  }
}
pracma::toc()

################################################################################
#################### Post-processing / summaries

results <- as.data.frame(results) %>% dplyr::distinct()

set_level <- results %>%
  dplyr::group_by(X1, X2, X3, topo_id, highly_robust) %>%
  dplyr::summarise(
    n_edges      = dplyr::n(),
    max_path_len = max(path_len, na.rm = TRUE),
    tot_intermed = sum(n_intermediate),
    all_paths    = paste(paste0(from_name, "->", to_name, ":", path),
                         collapse = " | "),
    .groups = "drop")

save(set_level, file = "set_level_SRC.RData")
write.csv(set_level, "set_level_SRC.csv", row.names = FALSE)
write.csv(results,   "edge_level_SRC.csv", row.names = FALSE)

high_rob_proteins <- unique(results[results$topo_id >= 35, c("X3", "topo_id")])
Unique_node_3     <- unique(results[, c("X3", "topo_id")])
not_high_rob_proteins <- results %>%
  dplyr::filter(!(X3 %in% high_rob_proteins$X3)) %>%
  dplyr::pull(X3) %>% unique()

cat("Unique rebound X3 matched (any topology): ",
    length(unique(Unique_node_3$X3)), "\n")
cat("Unique X3 in highly robust topologies:    ",
    length(unique(high_rob_proteins$X3)), "\n")
cat("Unique X3 only in non-highly-robust:       ",
    length(unique(not_high_rob_proteins)), "\n")

saveRDS(list(results = results, set_level = set_level,
             protein_name = protein_name, protein_index = protein_index),
        "SRC_result_with_paths.rds")