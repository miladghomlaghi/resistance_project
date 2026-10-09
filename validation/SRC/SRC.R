# rm(list = ls())
# gc()
library(dplyr)
load("../../protein_interactions_20_length_27_10_2022.RData")
load("../../Data/initial_info.RData")


## FOR SRC inhibitor



#################################################################################
#################################################################################
## finding the intersection between proteins exist in experimental data and and our network


my_data <- read.table("./pY221121.csv")

protein_index = NULL
protein_name_total = NULL

for (i in 1:dim(my_data)[1]) {
  #
  
  protein_name_total[i] <- strsplit(my_data[i, 1], '_')[[1]][2]

}
proteins_total_Exp_in_network <-
  unique(protein_name_total[protein_name_total %in%  nodes])

#################################################################################
#################################################################################
# calculating the distance between target and each rebound protein in the
# experimental data

# importing results of rebound analysis of experimental data
rebound_data <-
  read.csv(file = "bounceback_src_shared_two.csv",
           header = TRUE)

protein_distance = data.frame(col1 = character(),
                              col2 = numeric(),
                              stringsAsFactors = FALSE)

protein_name <- NULL

for (i in 1:dim(rebound_data)[1]) {
  
  
  protein_name[i] <- strsplit(rebound_data[i, 1], '_')[[1]][2]
}

protein_name<- unique(protein_name)

target_index <- which(nodes == "SRC")
target_index <- which(all_proteins_permut_proteins [, 1] ==target_index)


for (i in 1:length(protein_name)) {
  distance = 0
  
  print(i)
  protein_index <- which(nodes == protein_name[i])
  final_index = intersect(target_index,which(all_proteins_permut_proteins[, 2] == protein_index))
  
  if (!length(final_index) == 0) {
    
    distance <-
      length(all_proteins_permut_proteins[final_index, 3][[1]][[1]])
    # browser()
    protein_distance[nrow(protein_distance) + 1, ] <-
      c(protein_name[i], distance)
    
  }
}


protein_distance <- unique(protein_distance)
# total unique proteins that show rebound after treatment
unique_rebound_proteins <- unique(protein_name)
Disance_from_target<-as.integer(protein_distance[,2])

png(file="figures/SRC_Disance_from_target.png",
    width=200, height=350)
hist(Disance_from_target)
dev.off()

###################################################################################
# intersection between rebounded proteins and model predictions

# importing results related to topology search of SRC


load("../../Results/length_7_phospho/SRC.RData")

proteins_rebound_Exp_in_prediction <-
  unique(protein_name[protein_name %in%  output[[1]][, 3]])
length(proteins_rebound_Exp_in_prediction)
# unique node x3

length(unique(output[[1]][, 3]))

node3_proteins <- output[[1]]
test<-NULL

for(i in 1:length(proteins_rebound_Exp_in_prediction)){
  
  test[i]<-sum(proteins_rebound_Exp_in_prediction[i]==node3_proteins[, 3])
  
}


protein_name_intersect_nodes <- protein_name[protein_name %in% nodes]

proteins_rebound_Exp_not_in_prediction <-
  unique(protein_name[!(protein_name %in%  output[[1]][, 3])])

test1<-NULL
for(i in 1:length(proteins_rebound_Exp_not_in_prediction)){
  
  test1[i]<-sum(proteins_rebound_Exp_not_in_prediction[i]==node3_proteins[, 3])
  
}

node3_proteins[proteins_rebound_Exp_not_in_prediction[i]==node3_proteins[, 3],] 

################################################################################
## finding the total number of proteins with distance less than 20




library(igraph)


near_net <-
  igraph::ego(net,
              5,
              nodes = "SRC",
              mode = "all",
              mindist = 0)

near_net <- induced_subgraph(net, unlist(near_net))

near_nodes <- as_ids(V(near_net))




all_interaction_signor <-
  read.csv("validation/SRC/network_links_phospho_binding.csv")
output_table <- all_interaction_signor

for (i in 1:dim(all_interaction_signor)[1]) {
  if (grepl("up", all_interaction_signor[i, 3])) {
    output_table[i, 3] <- 1
    
  }
  
  if (grepl("down", all_interaction_signor[i, 3])) {
    output_table[i, 3] <- -1
    
  }
  
  
}
read.csv("network_links_phospho_binding_binary.csv")
length(unique(c(output_table[, 1], output_table[, 2])))


## What portion of node x3 are bindings in exp data and simulation?

intersect_sim_database <-
  all_interaction_signor[all_interaction_signor[, 2] %in% output[[1]][, 3] , c(2, 4)]
intersect_exp_database <-
  all_interaction_signor[all_interaction_signor[, 2] %in% proteins_total_Exp_in_network , c(2, 4)]


unique_exp <- unique(proteins_total_Exp_in_network)

binding_proteins = NULL

for (i in 1:length(unique_exp)) {
  if (length("phosphorylation" == all_interaction_signor[unique_exp[i] == all_interaction_signor[, 2], 4]) !=
      0) {
    binding_proteins =  c(binding_proteins, unique_exp[i])
    
  }
  
}

unique_exp <- unique(output[[1]][, 3])

# proteins as node 3 that are only interacting with binding in the database and there is not phosphorylation for them.

binding_proteins = NULL

for (i in 1:length(unique_exp)) {
  if (length("phosphorylation" == all_interaction_signor[unique_exp[i] == all_interaction_signor[, 2], 4]) ==
      0) {
    binding_proteins =  c(binding_proteins, unique_exp[i])
    
  }
  
}





##################################################################################
## finding interaction type of 3->1 and 3->2 for protein show rebound in the
## experimental data
library(data.table)

load("../../Results/length_20_revised/SRC.RData")

test <- output[[1]]
rebound_list <- unique(test[test[, 3] %in% protein_name, ])


second_col_proteins <- unique(rebound_list[, 2])

third_col_proteins <- unique(rebound_list[, 3])
node_index_one <- which(nodes == rebound_list[1, 1])

one_index = which(all_proteins_permut_proteins [, 2] ==
                    node_index_one)

node_index_two <- NULL
node_index_three <- NULL

for (j in 1:length(second_col_proteins)) {
  node_index_two[j] <- which(nodes == second_col_proteins[j])
  
}
for (j in 1:length(third_col_proteins)) {
  node_index_three[j] <- which(nodes == third_col_proteins[j])
  
}

# all_proteins_permut_values_t <- as.data.table(all_proteins_permut_values)
# all_proteins_permut_proteins_t <- as.data.table(all_proteins_permut_proteins)
# all_proteins_permut_values_t <- all_proteins_permut_values_t %>% mutate(V1 = as.numeric(V1),
#                                         V2 = as.numeric(V2))
# Set key columns for faster searching
# setkey(all_proteins_permut_values_t, V1, V2)
# 
# # Define your specific strings
# specific_string1 <- "your_string1"
# specific_string2 <- "your_string2"

# Filter the data.table based on specific strings in the first two columns
# result <- dt[column1 == specific_string1 & column2 == specific_string2, column3]

# tmp_sec_col_proteins_path <-
#   all_proteins_permut_values[intersect(
#     which(all_proteins_permut_proteins[, 2] %in% node_index_two),
#     which(all_proteins_permut_proteins[, 1] %in% node_index_three)
#   ), c(1, 4,5)]


tmp_thrd_col_proteins_path <-
  all_proteins_permut_values[intersect(one_index,
                                       which(all_proteins_permut_proteins[, 1] %in% node_index_three)), c(1,2, 3)]

tmp_thrd_col_proteins_first_path <-
  sapply(tmp_thrd_col_proteins_path[, 3], "[[", 1)


tmp_sec_col_proteins_path <-
  all_proteins_permut_values[intersect(one_index,
                                       which(all_proteins_permut_proteins[, 1] %in% node_index_two)), c(1,2, 3)]

tmp_sec_col_proteins_first_path <-
  sapply(tmp_sec_col_proteins_path[, 3], "[[", 1)

# three to Two node interaction
rebound_list_indices <- sapply(rebound_list[,c(2,3)], function(col) match(col, nodes))

tmp_sec_thrd_col_proteins_path <-
  all_proteins_permut_values[intersect(
    which(all_proteins_permut_values[, 2] %in% rebound_list_indices[,1]),
    which(all_proteins_permut_values[, 1] %in% rebound_list_indices[,2]))
  , c(1, 2,3)]

tmp_tmp <- as.data.frame(tmp_sec_thrd_col_proteins_path[,c(1,2)])
matching_indices <- which(apply(tmp_tmp, 1, function(row_B) {
  any(apply(as.data.frame(rebound_list_indices[c(1,2),]), 1, function(row_A) all.equal(row_A, row_B) == TRUE))
}))

rownames(rebound_list_indices) <- NULL
colnames(rebound_list_indices) <- NULL
rownames(tmp_tmp) <- NULL
colnames(tmp_tmp) <- NULL
rownames(all_proteins_permut_values) <- NULL
colnames(all_proteins_permut_values) <- NULL

# finding rows in all_proteins_permut_values matching with node 2 and 3 proteins 
# in  rebound_list_indices
a <- apply(rebound_list_indices, 1, function(row, row_index) {
  print(paste("Processing row index:", row_index))
  which((tmp_tmp[, 2] == row[1]) & (tmp_tmp[, 1] == row[2]))
}, row_index = 1:nrow(rebound_list_indices))

# finding shortest interaction type among the node 2 and 3
tmp_sec_thrd_col_proteins_first_path <-
  sapply(tmp_sec_thrd_col_proteins_path[a, 3], "[[", 1)

# finding how many node 3 protein do not interact with node 2. 
length(unique(rebound_list_indices[tmp_sec_thrd_col_proteins_first_path==0,2]))
length(intersect((tmp_thrd_col_proteins_path[,1]),(unique(rebound_list_indices[tmp_sec_thrd_col_proteins_first_path==0,2]))))

###############################################################################
# in the found proteins sets check between node three and two. extract ones with
# no link

rebound_list_mat <- as.matrix(rebound_list_indices)
all_proteins_permut_proteins_mat <- (all_proteins_permut_proteins[,c(1,2)])
# matching_indices <- which(rowSums(rebound_list_mat == all_proteins_permut_proteins_mat) == ncol(rebound_list_mat))


library(data.table)

# Convert matrices to data.tables
setDT(all_proteins_permut_proteins_mat)
setDT(as.data.frame(rebound_list_indices))

# Find matching indices
matching_indices <- all_proteins_permut_proteins_mat[, .I[rowSums(.SD == as.list(rebound_list_indices)) == ncol(rebound_list_indices)], by = 1:nrow(all_proteins_permut_proteins_mat)]$V1

# Display the indices
print(matching_indices)



################################################################################
# plotting distribution of output nodes from target in simulated
load("Results/length_20_revised/SRC.RData")

protein_distance = data.frame(
  col1 = character(), col2 = numeric(),stringsAsFactors = FALSE)
protein_name1=NULL

target_index <- which(nodes == "SRC")
target_index <- which(all_proteins_permut_proteins [, 1] ==target_index)


output_proteins<-unique(test[,3])

result <-foreach(i = 1:length(output_proteins), .combine = 'rbind') %dopar%{
  
  distance=0
  
  protein_index <- which(nodes == output_proteins)
  final_index = intersect(target_index, which(all_proteins_permut_proteins[, 2] == protein_index))
  
  if (!length(final_index)==0) {
    distance<- length(all_proteins_permut_proteins[final_index, 3][[1]][[1]])
    # browser()
    c(output_proteins,distance)
  }
  
}

# png(file="../../figures/SRC_Disance_from_target_simulations.png",
    # width=200, height=350)
# hist(Disance_from_target)
dev.off()

