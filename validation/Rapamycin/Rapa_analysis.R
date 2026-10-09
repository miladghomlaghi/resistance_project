# rm(list = ls())
# gc()
library(parallel)
library(foreach)
load("protein_interactions_20_length_27_10_2022.RData")
load("Data/initial_info.RData")

#################################################################################
#################################################################################
## Intersection between phosphosites measured and network proteins

all_Rapamycin_phospho_Data <-
  read.csv(file = "Rapamycin.csv",
           header = TRUE)

total_phosphosite_name<-(all_Rapamycin_phospho_Data[,2])
sum(total_phosphosite_name %in% nodes ) 


#################################################################################
#################################################################################
## finding distance between the nodes in the experimental data and target
exp_data_rebound <-
  read.csv(file = "Rapamycin_karina_constraint.csv",
           header = TRUE)

exp_data_rebound_in_network<-exp_data_rebound[exp_data_rebound[,1] %in%  nodes,1]

protein_distance = data.frame(col1 = character(),
                              col2 = numeric(),
                              stringsAsFactors = FALSE)

length(unique(exp_data_rebound_in_network[exp_data_rebound_in_network %in%  output[[1]][, 3]]))
length(unique(exp_data_rebound[exp_data_rebound[,1] %in%  output[[1]][, 3],]))

target_index <- which(nodes == "MTOR")
target_index <- which(all_proteins_permut_proteins [, 1] ==target_index)

for (i in 1:dim(my_data)[1]) { 
  
  print(i)
  protein_name1[i] <- my_data[i, 1]
}

protein_name1 <-unique(protein_name1)
for (i in 1:length(protein_name1)) { #
  distance=0
  
  protein_index <- which(nodes == protein_name1[i])
  final_index = intersect(target_index, which(all_proteins_permut_proteins[, 2] == protein_index))
  
  if (!length(final_index)==0) {
    distance<- length(all_proteins_permut_proteins[final_index, 3][[1]][[1]])
    # browser()
    protein_distance[nrow(protein_distance)+1,] <- c(protein_name1[i],distance)
  }
  
}


protein_distance <- unique(protein_distance)
# total unique proteins that show rebound after treatment
Disance_from_target<-as.integer(protein_distance[,2])

png(file="figures/MTOR_Disance_from_target.png",
    width=200, height=350)
hist(Disance_from_target)
dev.off()


#phopshosite that we have in our network
protein_name <- protein_name1
sum(protein_name1 %in% nodes ) 

protein_name<-unique((unlist(protein_name)))
sum(protein_name %in% nodes ) 



##################################################################################
## finding interaction type of 3->1 and 3->2 for protein show rebound in the
## experimental data

load("../../Results/length_7_phospho/MTOR.RData")

rebound_data <- output[[1]]
# rebound_data <- unique(test[test[,3] %in% protein_name,])

#####
second_col_proteins <- unique(rebound_data[, 2])

third_col_proteins <- unique(rebound_data[, 3])
node_index_one <- which(nodes == rebound_data[1, 1])

one_index = which(all_proteins_permut_proteins [, 2] ==
                    node_index_one)
###########
# finding index of protein in the network
node_index_two<-NULL
node_index_three<-NULL


for (j in 1:length(second_col_proteins)) {
  node_index_two[j] <- which(nodes == second_col_proteins[j])
  
}
for (j in 1:length(third_col_proteins)) {
  node_index_three[j] <- which(nodes == third_col_proteins[j])
  
}


## finding pathes from proteins in column three to two
tmp_sec_col_proteins_path <-
  all_proteins_permut_values[intersect(which(all_proteins_permut_proteins[, 2] %in% node_index_two),
                                       which(all_proteins_permut_proteins[, 1] %in% node_index_three)), c(1,4,5)]

## check if any unique protein in node 3 has no interaction with potential node 2 proteins
nodes3_interaction_type <- unlist(tmp_sec_col_proteins_path[,2])| unlist(tmp_sec_col_proteins_path[,3])
nodes3_interaction_table  <- cbind(unlist(tmp_sec_col_proteins_path[,1]) , as.integer(nodes3_interaction_type==FALSE))
colnames(nodes3_interaction_table)<-c("node","values")
nodes3_interaction_table<-data.frame(nodes3_interaction_table)

by_cyl  <- nodes3_interaction_table  %>% group_by(node)
by_cyl1 <- by_cyl %>% summarise(values = sum(values))



## finding paths from proteins in column three to one
tmp_thrd_col_proteins_path <-
  all_proteins_permut_values[intersect(one_index,
                                       which(all_proteins_permut_proteins[, 1] %in% node_index_three)), c(1,3)]

tmp_thrd_col_proteins_first_path <-
  sapply(tmp_thrd_col_proteins_path[,2], "[[", 1)

# Finding examples that are exactly similar to the core network motif

nodes3_interaction_table_three_one<-cbind(unlist(tmp_thrd_col_proteins_path[,1]),cbind(tmp_thrd_col_proteins_first_path,by_cyl1[,2]))
final<-nodes3_interaction_table_three_one[nodes3_interaction_table_three_one[,2]!=1 & nodes3_interaction_table_three_one[,3]!=0,]


length(unique(final[,1]))


###################################################################################


three_to_two <-
foreach(i = 1:100, .combine = rbind) %dopar% {# length(third_col_proteins)
    # for (i in 1:1){
    
    node_index_two <- NULL
    node_index_three <- NULL
    
    protein_col_two_specifc_to_three <- rebound_data[rebound_data[,3]==third_col_proteins[i],2]
    
    for (j in 1:length(protein_col_two_specifc_to_three)) {
      node_index_two[j] <- which(nodes == protein_col_two_specifc_to_three[j])
      
    }
    node_index_three <- which(nodes == third_col_proteins[i])
    
    three_index = which(all_proteins_permut_proteins [, 1] ==
                        node_index_three)
    
    tmp_thrd_col_proteins_path <-
      all_proteins_permut_values[intersect(three_index,
                                           which(all_proteins_permut_proteins[, 2] %in% node_index_two)), c(1,2,4,5)]
    
    tmp_thrd_col_proteins_path[!(unlist(tmp_thrd_col_proteins_path[,3]) | unlist(tmp_thrd_col_proteins_path[,4])),c(1,2)]
    
  }


restrict_MTOR<-three_to_two
#######################################################################
## finding paths from proteins in column three to two

test<- cbind(match(rebound_data[,3], nodes),match(rebound_data[,2], nodes))


result <- match(paste(test[,1], test[,2]), paste(all_proteins_permut_proteins[,1], all_proteins_permut_proteins[,2]))

three_to_two_semi_final = all_proteins_permut_values[result,c(1,2,4,5)]

three_to_two_index <- which(!(unlist(three_to_two_semi_final[,3]) | unlist(three_to_two_semi_final[,4])))

#######################################################
## finding paths from proteins in column three to one

one_index = which(all_proteins_permut_proteins [, 2] ==
                    node_index_one)


test<- cbind(match(rebound_data[,3], nodes),match(rebound_data[,1], nodes))


result_three_two_one <- match(paste(test[,1], test[,2]), paste(all_proteins_permut_proteins[,1], all_proteins_permut_proteins[,2]))


tmp_thrd_col_proteins_path <-
  all_proteins_permut_values[intersect(one_index,
                                       which(all_proteins_permut_proteins[, 1] %in% node_index_three)), c(1,3)]


tmp_thrd_col_proteins_first_path <-
  sapply(all_proteins_permut_values[result_three_two_one,3], "[[", 1)


index_three_two_one <- which(!(tmp_thrd_col_proteins_first_path==1))



new_rebound <- rebound_data[intersect(unlist(three_to_two_index),unlist(index_three_two_one)),]


################################################################################
# plotting distribution of output nodes from target in simulated
load("Results/length_20_revised/MTOR.RData")

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

png(file="figures/MTOR_Disance_from_target_simulations.png",
    width=200, height=350)
hist(Disance_from_target)
dev.off()



