# rm(list = ls())
# gc()

# load("protein_interactions_20_length_27_10_2022.RData")
# load("Data/initial_info.RData")
################################################################################
################################################################################
## Intersection between proteins measured and network proteins

all_MEK_phospho_Data <-
  read.csv(file = "./validation/MEK/AZD62.csv",
           header = TRUE)

total_protein_name <- (all_MEK_phospho_Data[, 2])
sum(total_protein_name %in% nodes)

################################################################################

## FOR MEK inhibitor

my_data <-
  read.csv(file = "./validation/SRC/AZD6244_karina_constraint.csv",
           header = TRUE)

protein_distance = data.frame(col1 = character(),
                              col2 = numeric(),
                              stringsAsFactors = FALSE)
protein_name1 = NULL


target_index <- which(nodes == "MAP2K1")
target_index <-
  which(all_proteins_permut_proteins [, 1] == target_index)

for (i in 1:dim(my_data)[1]) {
  print(i)
  protein_name1[i] <- my_data[i, 1]
}

protein_name1 <- unique(protein_name1)
for (i in 1:length(protein_name1)) {
  #
  distance = 0
  
  protein_index <- which(nodes == protein_name1[i])
  final_index = intersect(target_index,
                          which(all_proteins_permut_proteins[, 2] == protein_index))
  
  if (!length(final_index) == 0) {
    distance <-
      length(all_proteins_permut_proteins[final_index, 3][[1]][[1]])
    # browser()
    protein_distance[nrow(protein_distance) + 1, ] <-
      c(protein_name1[i], distance)
  }
  
}

# phosphosites show rebound in our network and exist in our network
sum(protein_name1 %in% nodes)
aa = protein_name1[protein_name1 %in% nodes]
# proteins show rebound in our network and exist in our network
length(unique(protein_name1[protein_name1 %in% nodes]))





protein_distance <- unique(protein_distance)

# total unique proteins that show rebound after treatment
Disance_from_target <- as.integer(protein_distance[, 2])

png(file = "figures/MEK_Disance_from_target.png",
    width = 200,
    height = 350)
hist(Disance_from_target)
dev.off()

# proteins that are in our network
protein_name <- unique(protein_name1)

##################################################################################
# distance of proteins that are recognized as rebounding proteins by Netscan

load("Results/length_20_revised/MAP2K1.RData")
protein_distance[protein_distance[,1] %in% output[[1]][,3],]
length(unique(output[[1]][,3]))

###############################################################################
# in network, rebound, unique
rebound_in_network_unique <- protein_name1[protein_name1 %in% nodes]
# in network, rebound, unique intersect examples

rebound_in_network_unique_intersect_examples <- (unique(rebound_in_network_unique[rebound_in_network_unique%in% output[[1]][,3]]))
protein_distance_rebound_in_network_unique_intersect_examples <- protein_distance[protein_distance[,1] %in% rebound_in_network_unique_intersect_examples,]

png(file = "figures/protein_distance_rebound_in_network_unique_intersect_examples.png",
    width = 200,
    height = 350)
hist(as.integer(protein_distance_rebound_in_network_unique_intersect_examples[,2]))
dev.off()

####################################
## finding interaction type of 3->1 and 3->2 for protein show rebound in the
## experimental data

test <- output[[1]]
rebound_data <- unique(test[test[, 3] %in% protein_name, ])

print("number of exp_proteins in network", sum(nodes %in% protein_name))
rebound_list <- rebound_data
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


tmp_sec_col_proteins_path <-
  all_proteins_permut_values[intersect(
    which(all_proteins_permut_proteins[, 2] %in% node_index_two),
    which(all_proteins_permut_proteins[, 1] %in% node_index_three)
  ), c(1, 4, 5)]


nodes3_interaction_type <-
  unlist(tmp_sec_col_proteins_path[, 2]) |
  unlist(tmp_sec_col_proteins_path[, 3])
nodes3_interaction_table  <-
  cbind(unlist(tmp_sec_col_proteins_path[, 1]) ,
        as.integer(nodes3_interaction_type == FALSE))
colnames(nodes3_interaction_table) <- c("node", "values")
nodes3_interaction_table <- data.frame(nodes3_interaction_table)

by_cyl  <- nodes3_interaction_table  %>% group_by(node)
by_cyl1 <- by_cyl %>% summarise(values = sum(values))


tmp_thrd_col_proteins_path <-
  all_proteins_permut_values[intersect(one_index,
                                       which(all_proteins_permut_proteins[, 1] %in% node_index_three)), c(1, 3)]

tmp_thrd_col_proteins_first_path <-
  sapply(tmp_thrd_col_proteins_path[, 2], "[[", 1)

nodes3_interaction_table_three_one <-
  cbind(
    unlist(tmp_thrd_col_proteins_path[, 1]),
    cbind(tmp_thrd_col_proteins_first_path, by_cyl1[, 2])
  )
final <-
  nodes3_interaction_table_three_one[nodes3_interaction_table_three_one[, 2] !=
                                       1 & nodes3_interaction_table_three_one[, 3] != 0, ]




###################################################################################

protein_name1 = NULL

target_index <- which(nodes == "MAP2K1")
target_index <-
  which(all_proteins_permut_proteins [, 1] == target_index)


output_proteins <- unique(test[, 3])

result <-
  foreach(i = 1:length(output_proteins), .combine = 'rbind') %dopar% { #
    # for (i in 1:100){
    
    
    distance = 0
    
    protein_index <- which(nodes == output_proteins[i])
    final_index = intersect(target_index,
                            which(all_proteins_permut_proteins[, 2] == protein_index))
    
    if (!length(final_index) == 0) {
      distance <-
        length(all_proteins_permut_proteins[final_index, 3][[1]][[1]])
      # browser()
      c(output_proteins[i], distance)
    }
    
  }

png(file = "figures/MEK_Disance_from_target_simulations.png",
    width = 200,
    height = 350)
hist(as.integer(result[,2]))
dev.off()
