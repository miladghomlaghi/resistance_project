# rm(list = ls())
# gc()
library(dplyr)
# load("../../protein_interactions_20_length_27_10_2022.RData")
load("../../Data/initial_info.RData")

load("new_initial_info_9_6_phospho_bind.RData")
load("../../protein_interactions_phospho_bind_7_6_2024.RData")


#################################################################################
#################################################################################
## raw Signor database

raw_signor_data <- read.csv("./7-6-24-phspho-binding.csv")

raw_signor_data_phospho_bind <- raw_signor_data%>%filter(MECHANISM%in% c("phosphorylation","dephosphorylation"))

data1 <- raw_signor_data_phospho_bind[, c("ENTITYA", "IDA")]
colnames(data1) <- c("ENTITY", "ID")

data2 <- raw_signor_data_phospho_bind[, c("ENTITYB", "IDB","RESIDUE")]
colnames(data2) <- c("ENTITY", "ID","RESIDUE")
#replacing Ser with S
data2 <- data2 %>%
  mutate(
    RESIDUE = gsub("(?i)Ser", "S", RESIDUE, perl = TRUE),
    RESIDUE = gsub("(?i)Thr", "T", RESIDUE, perl = TRUE),
    RESIDUE = gsub("(?i)Tyr", "Y", RESIDUE, perl = TRUE),
    tmp_pro_site = paste(ID) 
  )

#without considering phosphosites
UNIPORT_CODE_nodes <- distinct(rbind(data1,data2%>%select(c("ENTITY", "ID"))))

#with considering phosphosites
# UNIPORT_CODE_nodes <- distinct(rbind(data2))
# nodes <-unique(combined_data$ENTITY )

#################################################################################
#################################################################################
## finding the intersection between proteins exist in experimental data and and our network


my_data <- read.csv("./bounceback_src_total_valid.csv")
tmp_protein_code <- read.csv("./results.csv")%>%select(c("Phosphosite","Protein.ID","Residue.Both"))
colnames(my_data)[colnames(my_data) == "ProteinID"] <- "Phosphosite"
tmp_protein_code$Residue.Both <- sub(";.*", "", tmp_protein_code$Residue.Both)

my_data <- my_data %>%left_join(tmp_protein_code)%>%mutate(tmp_pro_site = paste(Protein.ID))

#without considering phosphosites
proteins_total_Exp_in_network <-
  unique(my_data%>%filter(Protein.ID %in%  UNIPORT_CODE_nodes$ID))

#with considering phosphosites
# proteins_total_Exp_in_network <-
#   unique(my_data%>%filter(tmp_pro_site %in%  UNIPORT_CODE_nodes$tmp_pro_site))

protein_index = NULL
protein_name_total = NULL

for (i in 1:dim(my_data)[1]) {
  
  protein_name_total[i] <- strsplit(my_data[i, 1], '_')[[1]][1]
  
}
tmp_protein_code <- read.csv("./results.csv")

unique(protein_name_total)
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
colnames(rebound_data)[colnames(rebound_data) == "ProteinID"] <- "Phosphosite"

rebound_data <- rebound_data %>%left_join(tmp_protein_code)%>%mutate(tmp_pro_site = paste(Protein.ID,Residue.Both))


protein_distance = data.frame(col1 = character(),
                              col2 = numeric(),
                              stringsAsFactors = FALSE)

protein_name <- NULL

for (i in 1:dim(rebound_data)[1]) {
  
  
  protein_name[i] <- strsplit(rebound_data[i, 1], '_')[[1]][1]
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
# node3_proteins: all predicted node 3 proteins
# UNIPORT_CODE_nodes: all proteins in the signor database with uniport code
# rebound_data: all rebounded proteins in the experimental data (without considering if they are in the network)
# tmp_protein_code: protein name-UNIPORT-residue from the SRC dataset

load("../../Results/length_7_phospho_bind_7_6/SRC.RData")


node3_proteins <- unique(output[[1]][, 3])
node3_uniport_code <-
  unique(UNIPORT_CODE_nodes$ID[UNIPORT_CODE_nodes$ENTITY %in%  node3_proteins])

#rebounded data exist in network considering protein:
rebound_data_in_network1 <-
  unique(rebound_data%>%filter(Protein.ID %in%  UNIPORT_CODE_nodes$ID)) # 472

proteins_rebound_Exp_in_prediction <-
  unique(rebound_data_in_network1$Protein.ID[rebound_data_in_network1$Protein.ID %in%  node3_uniport_code])

total_rebound_protein_intersect_nodes <-
  unique(rebound_data$Protein.ID[rebound_data$Protein.ID %in%  unique(UNIPORT_CODE_nodes$ID)])

length(total_rebound_protein_intersect_nodes) # 212 proteins rebounded and exists in our signor network
length(proteins_rebound_Exp_in_prediction) # among 212 there are 123 proteins rebounded and predicted by netscan
length(unique(proteins_total_Exp_in_network[proteins_total_Exp_in_network %in%  node3_proteins])) # denominator of raw score

length(unique(rebound_data$Protein.ID))



# unique node x3

length(unique(output[[1]][, 3]))

node3_proteins <- output[[1]]
unique(node3_proteins[, 3])


test<-NULL

for(i in 1:length(proteins_rebound_Exp_in_prediction)){
  
  test[i]<-sum(proteins_rebound_Exp_in_prediction[i]==node3_proteins[, 3])
  
}


protein_name_intersect_nodes <- protein_name[protein_name %in% nodes]
total_protein_name_intersect_nodes <- proteins_total_Exp_in_network[proteins_total_Exp_in_network %in% nodes]
no_rebound_exp_nodes_intersect <- proteins_total_Exp_in_network[!(proteins_total_Exp_in_network %in% total_rebound_protein_intersect_nodes)]

no_rebound_predict <- unique(UNIPORT_CODE_nodes$ID)[!(unique(UNIPORT_CODE_nodes$ID) %in% node3_uniport_code)]
no_rebound_exp_nodes_intersect[no_rebound_exp_nodes_intersect %in% no_rebound_predict]

proteins_rebound_Exp_not_in_prediction <-
  unique(protein_name_intersect_nodes[!(protein_name_intersect_nodes %in%  output[[1]][, 3])])



test1<-NULL
for(i in 1:length(proteins_rebound_Exp_not_in_prediction)){
  
  test1[i]<-sum(proteins_rebound_Exp_not_in_prediction[i]==node3_proteins[, 3])
  
}

node3_proteins[proteins_rebound_Exp_not_in_prediction[i]==node3_proteins[, 3],] 


a = unique(output[[1]][, 3])
a[!(a %in%  nodes )]



###################################################################################
#relationship between number of X2 and occurance of rebound
###################################################################################


SRC_predictions = as.data.frame(output[[1]])
colnames(SRC_predictions) = c("One","Two","Three")
colnames(UNIPORT_CODE_nodes) = c("Three","ID","Three")


X3_X2num = SRC_predictions%>%group_by(Three)%>%mutate(num_X2 = n())%>%select("Three",num_X2)%>%distinct()%>%arrange(desc(num_X2))
X3_X2num = X3_X2num %>% left_join(UNIPORT_CODE_nodes %>% select("Three", "ID"),by="Three")


# proteins in the experimental data that do not show rebound and exist in the network
no_rebound_exp_nodes_intersect <- proteins_total_Exp_in_network[!(proteins_total_Exp_in_network %in% rebound_data_in_network1)]%>%pull("Protein.ID")%>%unique()

#X2 number for non rebound
X2_num_non_rebound = X3_X2num%>%filter(ID %in% no_rebound_exp_nodes_intersect)%>%distinct()
# proteins in the experimental data that show rebound and exist in the network

rebound_exp_nodes_intersect = (rebound_data_in_network1%>%pull("Protein.ID")%>%unique())
X2_num_rebound = X3_X2num%>%filter(ID %in% rebound_exp_nodes_intersect)%>%distinct()

X2_num_non_rebound$Group <- "Non rebound"
X2_num_rebound$Group <- "Rebound"

combined <- bind_rows(X2_num_non_rebound, X2_num_rebound)

library(ggplot2)
library(ggpubr)

library(ggplot2)
library(ggpubr)

ggplot(combined, aes(x = Group, y = num_X2)) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    color = "grey20",
    fill = "grey92",
    linewidth = 0.7
  ) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.5, color = "grey30") +
  stat_compare_means(method = "t.test", label = "p.format", size = 5) +
  labs(
    x = NULL,
    y = "Number of X2 proteins",
    title = "Comparison of X2_num Between Protein Groups"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.line = element_line(color = "grey20", linewidth = 0.6),
    axis.ticks = element_line(color = "grey20", linewidth = 0.6),
    plot.margin = margin(10, 14, 10, 10)
  )

# comparing the number of X2 between these two groups


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
  read.csv("7-6-24-phspho-binding.csv")
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

load("../../Results/length_7_phospho_bind_7_6/SRC.RData")

test <- output[[1]]
rebound_list <- unique(test[test[, 3] %in% nodes, ])


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

