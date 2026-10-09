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
# saveRDS(c(list(targets),list(topos_to_match)),file="topos_to_match.rds")
tmp<-readRDS("topos_to_match.rds")
targets <- tmp[[1]]
topos_to_match <- tmp[[2]]


load("protein_interactions_15_4_2024.RData")

all_proteins_permut_proteins = all_proteins_permut_proteins[which(!(unlist(all_proteins_permut_values[,1])==unlist(all_proteins_permut_values[,2]))),]
all_proteins_permut_values = all_proteins_permut_values[which(!(unlist(all_proteins_permut_values[,1])==unlist(all_proteins_permut_values[,2]))),]

####################################################################################
############## load separate protein_interaction files and create variables ########

table_tmp <-   read.csv(file = "protein_interaction_phospho_4_24.csv",
                        header = TRUE)
nodes <- sort(unique(c(table_tmp[, 1], table_tmp[, 2])))
nodes <- unname(sort(nodes))
protein_range <- 1:length(nodes)
links<-table_tmp


protein1 = "FGFR1"
protein2 = "AKT2"
  
protein1_index=  which(nodes==protein1)
protein2_index=  which(nodes==protein2)

all_proteins_permut_proteins[(all_proteins_permut_proteins[,1] == protein1_index& all_proteins_permut_proteins[,2] == protein2_index),]









