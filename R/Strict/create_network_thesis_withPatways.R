#### NETWORK THESIS 05_MARCH

library(magrittr)
library(igraph)


functions_names<- list.files("./Strict/Data/signor_cancer", full.names = T )
pathways<-lapply(functions_names, read.delim)
pathways<-lapply(pathways, function(x){if(sum(is.element(c("PATHWAY_ID","PATHWAY_NAME"),names(x)))>1){
  x$PMID<-as.character(x$PMID)
  return(x[,-c(1,2)])
}else{
  x$PMID<-as.character(x$PMID)
  return(x)}})
  
pathways <-dplyr::bind_rows(pathways)

data_Signor<-read.delim("./Strict/Data/all_data_05_03_21.tsv")
data_Signor <-data_Signor[data_Signor$TAX_ID==9606,]


data_Signor<-rbind(data_Signor,pathways)

data_Targets <-unique(data_Signor[is.element(data_Signor$TYPEA, c("antibody","chemical")),"ENTITYB"])


data_Signor <-data_Signor[data_Signor$EFFECT != "unknown",]
data_Signor <-data_Signor[data_Signor$EFFECT != "form complex",]
data_Signor <-data_Signor[!is.na(data_Signor$EFFECT),]
data_Signor <-data_Signor[!is.na(data_Signor$SCORE),]


sum(is.na(data_Signor$SCORE))
data_Signor[5848,]

levels(factor(data_Signor$TYPEA))
levels(factor(data_Signor$MECHANISM))


data_Signor <-data_Signor[is.element(data_Signor$TYPEA, c("complex","fusion protein", "protein", "smallmolecule")),]
data_Signor <-data_Signor[is.element(data_Signor$TYPEB, c("complex","fusion protein", "protein", "smallmolecule")),]
#data_Signor <-data_Signor[is.element(data_Signor$MECHANISM, c("phosphorylation","dephosphorylation")),]



data_Signor <-data_Signor[data_Signor$DIRECT=="YES",]


levels(factor(data_Signor$MECHANISM))



data_Signor$EFFECT_BIN<-ifelse(grepl(c("down"),data_Signor$EFFECT),-1,1)

data_Signor<- data_Signor[,c("ENTITYA","ENTITYB","EFFECT_BIN","EFFECT",
                                     colnames(data_Signor)[!is.element(colnames(data_Signor),
                                                                           c("ENTITYA","ENTITYB","EFFECT_BIN", "EFFECT"))])]

data_Signor<-data_Signor%>% dplyr::distinct(ENTITYA, ENTITYB, .keep_all = TRUE)


tmp_nodes<-unique(c(data_Signor$ENTITYA,data_Signor$ENTITYB))
tmp_links<-data_Signor
tmp_net <- igraph::graph_from_data_frame(tmp_links,vertices = tmp_nodes,directed = T)



descompose_net<- igraph::decompose(tmp_net)
connected_net <- descompose_net[[which((igraph::components(tmp_net)$csize==max(igraph::components(tmp_net)$csize)))]]

nodes <- as_ids(V(connected_net))
net <- connected_net
links <- tmp_links
targets <- data_Targets[is.element(data_Targets,nodes)]


save(links,nodes,net,file = './Strict/Data/network_thesis_withPathways.Rdata')




#data_String<-read.delim("./Strict/Data/9606.protein.links.full_05042021.txt", sep = ' ')
