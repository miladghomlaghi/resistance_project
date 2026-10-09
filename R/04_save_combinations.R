nodes_num <- 3
topology <-8213
network_info <- "signor_PathwaysPlus" 

load(paste0("./Data/matching/",network_info,".RData", collapse = ''))
load("./Data/matching/3Nodes/signor_PathwaysPlus/merged_signor_PathwaysPlus.RData")
igraph::V(net)$name[igraph::V(net)$name =="1-phosphatidyl-1D-myo-inositol 3,4,5-trisphosphate"]<- "PIP3"
links[links=="1-phosphatidyl-1D-myo-inositol 3,4,5-trisphosphate"]<- "PIP3"
nodes[nodes=="1-phosphatidyl-1D-myo-inositol 3,4,5-trisphosphate"]<- "PIP3"

combinations<- get_combinations(merged_data,	topology, nodes_num)
#test <-frequency_match[126,3]
combinations[combinations=="1-phosphatidyl-1D-myo-inositol 3,4,5-trisphosphate"]<- "PIP3"

net_sizes<-get_subnetworks_size(combinations,net, nodes_num)
combinations <- net_sizes[order(net_sizes$NumNodes),]


combinations<-cbind("Combi_ID"=paste0(paste0("T",topology,"."),seq(1,nrow(combinations))),combinations)
save(file=paste0("./Data/matching/3Nodes/",network_info,"/T",topology,"_Networks.RData"),combinations)
write_xlsx(combinations, paste0("./Data/matching/3Nodes/",network_info,"/T",topology,"_Networks.xlsx"))

#save(file="testFrequency.RData",frequency_match )

# 
# T6431<-cbind("Combi_ID"=paste0("T6431.",seq(1,nrow(T6431))),T6431)
# T6434<-cbind("Combi_ID"=paste0("T6434.",seq(1,nrow(T6434))),T6434)
# 
# write.xlsx(T40847, file="combinationsSignor.xlsx", sheetName = "T40847", 
#            col.names = TRUE, row.names = FALSE, append = FALSE)
# write.xlsx(T5837, file="combinationsSignor.xlsx", sheetName = "T5837", 
#            col.names = TRUE, row.names = FALSE, append = TRUE)
# write.xlsx(T6434, file="combinationsSignor.xlsx", sheetName = "T6434", 
#            col.names = TRUE, row.names = FALSE, append = TRUE)
# 
# 
# test3<-T5837
# test_dis <- test3%>%distinct(NodeA, NodeC, .keep_all = TRUE)
# test4 <- test_dis[test_dis$NodeA %in% targetsSignor,]
# 
# write.xlsx(test4, file="combinationsSignor_filter.xlsx", sheetName = "T8213", 
#            col.names = TRUE, row.names = FALSE, append = FALSE)
# write.xlsx(test4, file="combinationsSignor_filter.xlsx", sheetName = "T5837", 
#            col.names = TRUE, row.names = FALSE, append = TRUE)
# 
# 
# nodes_combi<- merged_data[merged_data[,1]==28235,3]
# tmp <- nodes_combi
# tmp_split <- unlist(strsplit(tmp,"/"))
# #grep("EGFR",tmp_split)
# tmp_split<-tmp_split[grep("EGFR",tmp_split)]
# tmp_split<-tmp_split[grep("STAT3",tmp_split)]

