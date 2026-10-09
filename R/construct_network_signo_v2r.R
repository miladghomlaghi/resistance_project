### creating networks from signor database
library(magrittr)
library(igraph)
files_folder <- "Cancer_Pathways"

matching_Datafiles <- list.files(path=paste0("./RowData/Databases/Signor/",files_folder,collapse =''), 
                                 full.names = T)

tmp_pathways <- dplyr::bind_rows(lapply(matching_Datafiles, read.delim))

tmp_pathways<- data.frame(lapply(tmp_pathways, function (x) {
  gsub("1-phosphatidyl-1D-myo-inositol 3,4,5-trisphosphate","PI3P",x)}))

tmp_pathways<- data.frame(lapply(tmp_pathways, function (x) {
  gsub("3',5'-cyclic AMP","cAMP",x)}))


signor_pathways <-tmp_pathways
signor_pathways <-signor_pathways[signor_pathways$EFFECT != "unknown",]
signor_pathways <-signor_pathways[signor_pathways$TYPEA != "phenotype",]
signor_pathways <-signor_pathways[signor_pathways$TYPEB != "phenotype",]
signor_pathways <- signor_pathways[signor_pathways$TAX_ID==9606,]
signor_pathways<-signor_pathways[!is.na(signor_pathways$TAX_ID),]
row.names(signor_pathways)<-1:nrow(signor_pathways)
signor_pathways$STATE <- "original"
#issue_proteins <- c("SMAD2","SMAD3","SMAD4","BAK","BAX","RBPJ","NOTCH","IKK","IKBKG","STAT1","STAT3","JAK")
issue_proteins <- "SMAD2|SMAD3|SMAD4|BAK|BAX|RBPJ|NOTCH|IKK|IKBKG|STAT1|STAT3|JAK"
data_complexes <- signor_pathways[unique(c(grep(issue_proteins,signor_pathways$ENTITYA),
                                           grep(issue_proteins,signor_pathways$ENTITYB))),]
no_issue_proteins <- signor_pathways[!(row.names(signor_pathways) %in% row.names(data_complexes)),]
data_complexes<-data_complexes %>%  dplyr::distinct(ENTITYA, ENTITYB,EFFECT, .keep_all = TRUE)
no_issue_proteins<-no_issue_proteins %>%  dplyr::distinct(ENTITYA, ENTITYB,EFFECT, .keep_all = TRUE)

######levels(data_complexes$EFFECT) <- c(levels(data_complexes$EFFECT), "up-complex")
######data_complexes$EFFECT[data_complexes$EFFECT=="form complex"]<-"up-complex"

data_complexes$EFFECT_BIN<-ifelse(grepl(c("down"),data_complexes[,'EFFECT']),-1,1)
no_issue_proteins$EFFECT_BIN<-ifelse(grepl(c("down"),no_issue_proteins[,'EFFECT']),-1,1)

data_complexes<-data_complexes%>% dplyr::distinct(ENTITYA, ENTITYB,EFFECT_BIN, .keep_all = TRUE)
no_issue_proteins<-no_issue_proteins%>% dplyr::distinct(ENTITYA, ENTITYB,EFFECT_BIN, .keep_all = TRUE)


row.names(data_complexes)<-1:nrow(data_complexes)

data_complexes[c(33,34,35),"ENTITYB"]<-"JAK1/STAT1/STAT3"
data_complexes[c(33,34,35),"STATE"]<-"modified"   #SRC activates STAT1 and STAT3
row.names(data_complexes)<-1:nrow(data_complexes)


data_complexes <- data_complexes[-c(19,23,28,29), ] #to remove SMAD4 we simplified by considering SMAD2 upregulating SMAD2|SMAD4 complexes
row.names(data_complexes)<-1:nrow(data_complexes)
data_complexes[48,"ENTITYB"] <- "SMAD3/SMAD4" #SMAD6 acting directly in SMAD3/SMAD4
data_complexes[52,"ENTITYB"] <- "SMAD2/SMAD4"
data_complexes[c(48,52),"STATE"]<-"modified"
 

data_complexes[56,"ENTITYA"] <- "SMURF2" 
data_complexes[56,"STATE"]<-"modified"
data_complexes <- data_complexes[-53, ] 


signor_pathways<-rbind(data_complexes,no_issue_proteins)

signor_pathways<- signor_pathways[,c("ENTITYA","ENTITYB","EFFECT_BIN",colnames(signor_pathways)[!is.element(colnames(signor_pathways),c("ENTITYA","ENTITYB","EFFECT_BIN"))])]

modelling_network <- signor_pathways
modelling_network$SOURCE <- "signor"
row.names(modelling_network)<-1:nrow(modelling_network)
modelling_network<-modelling_network[-c(129, 130, 131,133),] #removing STK3/4, MOB1, Lats families
row.names(modelling_network)<-1:nrow(modelling_network)


modelling_network <-modelling_network[modelling_network$TYPEB != "stimulus",]
modelling_network <-modelling_network[modelling_network$TYPEA != "stimulus",]



nodes<-unique(c(modelling_network$ENTITYA,modelling_network$ENTITYB))
links<-modelling_network
net= graph_from_data_frame(links,vertices = nodes,directed = T)
save(links,nodes,net,file = 'cancer_network.Rdata')

which(nodes=='STAT3')


# levels(modelling_network$ENTITYB) <- c(levels(modelling_network$ENTITYB), "YAP1/TEAD")
# levels(modelling_network$MECHANISM) <- c(levels(modelling_network$MECHANISM), "unknown")
# modelling_network[186,"ENTITYB"]<-"YAP1/TEAD"
# modelling_network[186,"TYPEB"]<-"complex"
# modelling_network[186,"STATE"]<-"modified"
# modelling_network[186,"SOURCE"]<-"cellSignaling"
# modelling_network<-rbind(modelling_network,c('AKT','mTORC1',-1,'modified','down-regulates','unknown','proteinfamily','complex',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('AKT','TSC1/TSC2',-1,'modified','down-regulates','unknown','proteinfamily','complex',"cellSignaling"))
# levels(modelling_network$ENTITYB) <- c(levels(modelling_network$ENTITYB), "CRAF")
# modelling_network<-rbind(modelling_network,c('AKT','CRAF',1,'modified','up-regulates','unknown','proteinfamily','protein',"cellSignaling"))
# levels(modelling_network$ENTITYB) <- c(levels(modelling_network$ENTITYB), "PFK")
# modelling_network<-rbind(modelling_network,c('AKT','PFK',1,'modified','up-regulates','unknown','proteinfamily','proteinfamily',"cellSignaling"))
# levels(modelling_network$ENTITYB) <- c(levels(modelling_network$ENTITYB), "BAD")
# modelling_network<-rbind(modelling_network,c('14-3-3','BAD',-1,'modified','down-regulates','binding','proteinfamily','protein',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('AKT','XIAP',1,'modified','up-regulates','unknown','proteinfamily','protein',"cellSignaling"))
# 
# modelling_network<-rbind(modelling_network,c('AMPK','TSC1/TSC2',1,'modified','up-regulates','unknown','protein','complex',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('AMPK','PFK',1,'modified','up-regulates','unknown','protein','proteinfamily',"cellSignaling"))
# levels(modelling_network$ENTITYB) <- c(levels(modelling_network$ENTITYB), "BCLxL")
# levels(modelling_network$ENTITYA) <- c(levels(modelling_network$ENTITYA), "BAD")
# modelling_network<-rbind(modelling_network,c('BAD','BCLxL',-1,'modified','down-regulates','unknown','protein','protein',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('BAX','BCL2',-1,'modified','down-regulates','unknown','protein','protein',"cellSignaling"))
# 
# levels(modelling_network$ENTITYA) <- c(levels(modelling_network$ENTITYA), "BMPI/BMPII", "SMAD1/SMAD4","SMAD5/SMAD4","SMAD8/SMAD4",
#                                        "SMAD1","SMAD5","SMAD8")
# levels(modelling_network$ENTITYB) <- c(levels(modelling_network$ENTITYB), "SMAD1","SMAD5","SMAD8","BMPI/BMPII", 
#                                        "SMAD1/SMAD4","SMAD5/SMAD4","SMAD8/SMAD4")
# modelling_network<-rbind(modelling_network,c('BMPI/BMPII','SMAD1',1,'modified','up-regulates','unknown','complex','protein',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('BMPI/BMPII','SMAD5',1,'modified','up-regulates','unknown','complex','protein',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('BMPI/BMPII','SMAD8',1,'modified','up-regulates','unknown','complex','protein',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('BMPI/BMPII','SMAD8',1,'modified','up-regulates','unknown','complex','protein',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('SMAD1','SMAD1/SMAD4',1,'modified','up-regulates','unknown','protein','complex',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('SMAD5','SMAD5/SMAD4',1,'modified','up-regulates','unknown','protein','complex',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('SMAD8','SMAD8/SMAD4',1,'modified','up-regulates','unknown','protein','complex',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('SMAD6','SMAD1/SMAD4',-1,'modified','down-regulates','unknown','protein','complex',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('SMURF','SMAD1/SMAD4',-1,'modified','down-regulates','unknown','protein','complex',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('SMAD6','SMAD5/SMAD4',-1,'modified','down-regulates','unknown','protein','complex',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('SMURF','SMAD5/SMAD4',-1,'modified','down-regulates','unknown','protein','complex',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('SMAD6','SMAD8/SMAD4',-1,'modified','down-regulates','unknown','protein','complex',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('SMURF','SMAD8/SMAD4',-1,'modified','down-regulates','unknown','protein','complex',"cellSignaling"))
# 
# modelling_network<-rbind(modelling_network,c('SMAD1/SMAD4','SMAD7',1,'modified','up-regulates','unknown','complex','protein',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('SMAD5/SMAD4','SMAD7',1,'modified','up-regulates','unknown','complex','protein',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('SMAD8/SMAD4','SMAD7',1,'modified','up-regulates','unknown','complex','protein',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('SMAD7','BMPI/BMPII',-1,'modified','down-regulates','unknown','protein','complex',"cellSignaling"))
# 
# 
# levels(modelling_network$ENTITYA) <- c(levels(modelling_network$ENTITYA), "CAMK2A","RASGRP1","RASGRF1","SYNGAP1")
# levels(modelling_network$ENTITYB) <- c(levels(modelling_network$ENTITYB), "RASGRP1","RASGRF1","SYNGAP1")
# modelling_network<-rbind(modelling_network,c('CAMK2A','RASGRP1',1,'modified','up-regulates','unknown','protein','protein',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('CAMK2A','RASGRF1',1,'modified','up-regulates','unknown','protein','protein',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('CAMK2A','SYNGAP1',1,'modified','up-regulates','unknown','protein','protein',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('RASGRP1','HRAS',1,'modified','up-regulates','unknown','protein','protein',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('RASGRF1','HRAS',1,'modified','up-regulates','unknown','protein','protein',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('SYNGAP1','HRAS',-1,'modified','down-regulates','unknown','protein','protein',"cellSignaling"))
# 
# 
# levels(modelling_network$ENTITYB) <- c(levels(modelling_network$ENTITYB), "RAPGEF3")
# modelling_network<-rbind(modelling_network,c('3,5\'-cyclic AMP','RAPGEF3',1,'modified','up-regulates','unknown','protein','protein',"cellSignaling"))
# modelling_network<-rbind(modelling_network,c('CAMP','PRKACA',1,'modified','up-regulates','unknown','protein','protein',"cellSignaling"))


