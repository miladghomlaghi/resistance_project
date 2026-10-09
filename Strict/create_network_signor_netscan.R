#### SIGNOR FULL NETWORK NETSCAN

library(prodlim)
library(magrittr)
library(igraph)


data_Signor<-read.delim("./Strict/Data/all_data_05_03_21.tsv")
data_Signor <-data_Signor[data_Signor$TAX_ID==9606,]


data_Signor <-data_Signor[data_Signor$EFFECT != "unknown",]
data_Signor <-data_Signor[data_Signor$EFFECT != "form complex",]
data_Signor <-data_Signor[!is.na(data_Signor$EFFECT),]
data_Signor <-data_Signor[!is.na(data_Signor$SCORE),]


data_Signor <-data_Signor[is.element(data_Signor$TYPEA, c("complex","fusion protein", "protein", "smallmolecule")),]
data_Signor <-data_Signor[is.element(data_Signor$TYPEB, c("complex","fusion protein", "protein", "smallmolecule")),]


data_Signor <-data_Signor[data_Signor$DIRECT=="YES",]


data_Signor$EFFECT_BIN<-ifelse(grepl(c("down"),data_Signor$EFFECT),-1,1)

data_Signor<- data_Signor[,c("ENTITYA","ENTITYB","EFFECT_BIN","EFFECT",
                             colnames(data_Signor)[!is.element(colnames(data_Signor),
                                                               c("ENTITYA","ENTITYB","EFFECT_BIN", "EFFECT"))])]

data_Signor<-data_Signor%>% dplyr::distinct(ENTITYA, ENTITYB, .keep_all = TRUE)


data_Signor<-data_Signor[, c("ENTITYA","ENTITYB","EFFECT_BIN")]

data_Signor$Database<-"SIGNOR"

names(data_Signor)<-c("ENTITYA","ENTITYB","STATUS","DATABASE")

nodes<-
  unique(c(as.character(data_Signor$ENTITYA),as.character(data_Signor$ENTITYB)))
links<-data_Signor
net <- igraph::graph_from_data_frame(links,vertices = nodes,directed = T)


save(links,nodes,net,file = './Strict/Data/SIGNOR_NETSCAN.RData')


