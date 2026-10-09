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
load("protein_interactions_20_length_27_10_2022.RData")



result_table<-output[[1]]


potential_topo <- result_table[result_table[,3]=="AKT3",]


index_protein1<-which(nodes==potential_topo[1,1])
index_protein2<-which(nodes==potential_topo[1,2])
index_protein3<-which(nodes==potential_topo[1,3])


index_1_2 <- which(all_proteins_permut_proteins[,1]==index_protein1&
                     all_proteins_permut_proteins[,2]==index_protein2)

# Find the maximum length of the inner lists
max_length <- max(sapply(all_proteins_permut_proteins[index_1_2,3], length))

# Fill shorter lists with NA and convert to data frame
df <- t(as.data.frame(do.call(rbind, lapply(all_proteins_permut_proteins[index_1_2,3], function(x) {
  length(x) <- max_length
  return(x)
})), stringsAsFactors = FALSE))

df_1_2<- cbind(df,unlist(all_proteins_permut_values[index_1_2,3][1]))


index_1_3 <- which(all_proteins_permut_proteins[,1]==index_protein1&
                     all_proteins_permut_proteins[,2]==index_protein3)



# Find the maximum length of the inner lists
max_length <- max(sapply(all_proteins_permut_proteins[index_2_3,3], length))

# Fill shorter lists with NA and convert to data frame
df_1_3 <- t(as.data.frame(do.call(rbind, lapply(all_proteins_permut_proteins[index_1_3,3], function(x) {
  length(x) <- max_length
  return(x)
})), stringsAsFactors = FALSE))

df_1_3<- cbind(df_1_3,unlist(all_proteins_permut_values[index_1_3,3][1]))





index_2_3 <- which(all_proteins_permut_proteins[,1]==index_protein2&
                     all_proteins_permut_proteins[,2]==index_protein3)



# Find the maximum length of the inner lists
max_length <- max(sapply(all_proteins_permut_proteins[index_2_3,3], length))

# Fill shorter lists with NA and convert to data frame
df_2_3 <- t(as.data.frame(do.call(rbind, lapply(all_proteins_permut_proteins[index_2_3,3], function(x) {
  length(x) <- max_length
  return(x)
})), stringsAsFactors = FALSE))

df_2_3<- cbind(df_2_3,unlist(all_proteins_permut_values[index_2_3,3][1]))



