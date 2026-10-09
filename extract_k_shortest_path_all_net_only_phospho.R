closeAllConnections()
rm(list = ls())
gc()

library(dplyr)
library(foreach)
library(pracma)
library(igraph)
library(doParallel)

load("Data/initial_info.RData")
load("Data/permut.rData")
source("extract_k_shortest_function.R")
functions_names <- list.files("./data/functions", full.names = T)
lapply(functions_names, source)



###########################
table_tmp <-   read.csv(file = "protein_interaction_phospho_4_24.csv",
                        header = TRUE)
nodes <- sort(unique(c(table_tmp[, 1], table_tmp[, 2])))
nodes <- unname(sort(nodes))
protein_range <- 1:length(nodes)
dist <- 100

all_proteins_permut <-
  gtools::permutations(length(nodes),
                       2,repeats.allowed = F)
all_proteins_permut<-rbind(all_proteins_permut, cbind(c(1:length(nodes)),
                                                      c(1:length(nodes))))


links <- table_tmp
net <-
  igraph::graph_from_data_frame(table_tmp, directed = TRUE, vertices = NULL)

# browser()

steps <- seq(1, nrow(all_proteins_permut), 10000)
steps <- c(steps, nrow(all_proteins_permut))


numCores <- bigstatsr::nb_cores()

for (i in 14:length(steps)) {
  print(i)
  doParallel::registerDoParallel(numCores)
  
  
  tmp_proteins_permut <- all_proteins_permut[steps[i - 1]:steps[i],]
  
  path_proteins <- path_proteins_function(net,
                                          nodes,
                                          links,
                                          tmp_proteins_permut)
  
  save(
    path_proteins,
    file = paste0(
      "./Strict/Data/only_phospho/shortest_path_",
      i,
      ".RData",
      collapse = ""
    )
  )
  gc()
  doParallel::stopImplicitCluster()
  
  closeAllConnections()
  
}
