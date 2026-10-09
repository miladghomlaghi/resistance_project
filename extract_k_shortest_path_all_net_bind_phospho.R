closeAllConnections()
rm(list = ls())
gc()

library(dplyr)
library(foreach)
library(pracma)
library(igraph)
library(doParallel)

load("new_initial_info_9_6_phospho_bind.RData")
load("Data/permut.rData")
source("extract_k_shortest_function.R")
functions_names <- list.files("./data/functions", full.names = T)
lapply(functions_names, source)



###########################

protein_range <- 1:length(nodes)
dist <- 7

all_proteins_permut <-
  gtools::permutations(length(nodes),
                       2,repeats.allowed = F)
all_proteins_permut<-rbind(all_proteins_permut, cbind(c(1:length(nodes)),
                                                      c(1:length(nodes))))





##################################################################################
# main part
# browser()

steps <- seq(1, nrow(all_proteins_permut), 10000)
steps <- c(steps, nrow(all_proteins_permut))


numCores <- bigstatsr::nb_cores()

for (i in 2:length(steps)) {
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
      "./Strict/Data/phospho_bind_7_6_24/shortest_path_",
      i,
      ".RData",
      collapse = ""
    )
  )
  gc()
  doParallel::stopImplicitCluster()
  
  closeAllConnections()
  
}


routes <- list.files("./Strict/Data/phospho_bind_7_6_24/", full.names = T)

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
save(all_proteins_permut_proteins,all_proteins_permut_values,file = "protein_interactions_phospho_bind_7_6_2024.RData")






