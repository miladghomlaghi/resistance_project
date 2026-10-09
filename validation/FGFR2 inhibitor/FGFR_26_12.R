# rm(list = ls())
# gc()
library(dplyr)
# load("../../protein_interactions_20_length_27_10_2022.RData")
load("../../Data/initial_info.RData")

# load("new_initial_info_9_6_phospho_bind.RData")
load("../../protein_interactions_phospho_bind_7_6_2024.RData")


#################################################################################
#################################################################################
## raw Signor database

raw_signor_data <- read.csv("./7-6-24-phspho-binding.csv")

raw_signor_data_phospho_bind <- raw_signor_data%>%filter(MECHANISM%in% c("phosphorylation","dephosphorylation"))

data1 <- raw_signor_data_phospho_bind[, c("ENTITYA", "IDA")]
colnames(data1) <- c("ENTITY", "ID")

data2 <- raw_signor_data_phospho_bind[, c("ENTITYB", "IDB")]
colnames(data2) <- c("ENTITY", "ID")

UNIPORT_CODE_nodes <- distinct(rbind(data1, data2))
# nodes <-unique(combined_data$ENTITY )

#################################################################################
#################################################################################
## finding the intersection between proteins exist in experimental data and and our network


my_data <- read.csv("./bounceback_FGFR_total_valid.csv")
# tmp_protein_code <- read.csv("./CDK46_clean.csv")
# tmp_protein_code = tmp_protein_code[,c(1,2,3)]
names(my_data)[1:3] <- c("Protein.ID", "AminoAcid", "Phosphosite")


proteins_total_Exp_in_network <-
  unique(my_data$Protein.ID[my_data$Protein.ID %in%  UNIPORT_CODE_nodes$ID])


protein_index = NULL
protein_name_total = NULL

for (i in 1:dim(my_data)[1]) {
  
  protein_name_total[i] <- strsplit(my_data[i, 1], '_')[[1]][1]
  
}

unique(protein_name_total)
proteins_total_Exp_in_network <-
  unique(protein_name_total[protein_name_total %in%  nodes])


#################################################################################
################################################################################

###################################################################################
# intersection between rebounded proteins and model predictions

# importing results related to topology search of CDK

# total proteins in the signor database: 2681
# my_data: all proteins exist in the experimental data
# proteins_total_Exp_in_network: proteins in experimental data that exist in the 
# signor database: 635
# rebound_data: rebounded data: 2529 phosphosites 2529 proteins
# UNIPORT_CODE_nodes: UNIPORT code of all proteins in the signor 
# total_rebound_protein_intersect_nodes: total rebounded proteins that exist in signor database: 88
# node3_uniport_code: UNIPORT code of node 3 proteins 
# proteins_rebound_Exp_in_prediction: rebounded proteins that are in the network and 
# were predicted by Netscan

rebound_data <-
  read.csv(file = "bounceback_FGFR_shared_two.csv",
           header = TRUE)

rebound_data_in_network = rebound_data%>%filter(GeneNames%in% UNIPORT_CODE_nodes$GeneNames)

colnames(UNIPORT_CODE_nodes)=c("GeneNames","ID")
rebound_data_in_network = rebound_data_in_network%>%left_join(UNIPORT_CODE_nodes,by = "GeneNames")
colnames(rebound_data_in_network)[colnames(rebound_data_in_network) == "ID"] <- "Proteins"



load("../../Results/length_7_phospho_bind_7_6/FGFR2.RData")

total_rebound_protein_intersect_nodes <-
  unique(rebound_data_in_network$GeneNames[rebound_data_in_network$GeneNames %in%  unique(UNIPORT_CODE_nodes$GeneNames)])

node3_proteins <- unique(output[[1]][, 3])
node3_uniport_code <-
  unique(UNIPORT_CODE_nodes$ID[UNIPORT_CODE_nodes$GeneNames %in%  node3_proteins])

proteins_rebound_Exp_in_prediction <-
  unique(rebound_data_in_network$GeneNames[rebound_data_in_network$GeneNames %in%  node3_proteins])

proteins_rebound_Exp_not_in_prediction <-
  unique(rebound_data_in_network$GeneNames[!(rebound_data_in_network$GeneNames %in%  node3_proteins)])

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

