#### STRING FULL NETWORK NETSCAN


protein_action <- read.delim("RowData/Databases/String/9606.protein.actions.v11.0.txt")
protein.info <- read.delim("RowData/Databases/String/9606.protein.info.v11.0.txt", header=T, quote = "")

#protein.aliases <- read.delim("RowData/Databases/String/9606.protein.aliases.v11.0.txt", header=F, comment.char="#")
#colnames(protein.aliases)<-c("string_protein_id","alias","source")
#protein.aliases <- protein.aliases[grepl(c("kinase"), protein.aliases[,2]),]
#protein.links <- read.delim("RowData/Databases/String/9606.protein.links.full.v11.02.txt", header=T,sep=" ",row.names=NULL)

protein_action$ENTITYA <-""
protein_action$ENTITYB <-""
head(protein_action)
protein_action$ENTITYA<-protein.info[match(protein_action$item_id_a,protein.info$protein_external_id),2]
protein_action$ENTITYB<-protein.info[match(protein_action$item_id_b,protein.info$protein_external_id),2]
head(protein_action)

protein_action <- protein_action[protein_action$a_is_acting =="t",]
protein_action<-protein_action[protein_action$is_directional=="t",]
protein_action<-protein_action[is.element(protein_action$action, c("activation","inhibition")),]
protein_action<-protein_action[protein_action$score>900,]

###converting data to binary (1 activation, -1 inhibition).

protein_action$STATUS<- ifelse(protein_action$action=="activation",1,-1)

data_STRING <- protein_action[,c("ENTITYA","ENTITYB","STATUS")]
data_STRING$DATABASE<-"STRING"

nodes<-unique(c(as.character(data_STRING$ENTITYA),as.character(data_STRING$ENTITYB)))
links<-data_STRING
net <- igraph::graph_from_data_frame(links,vertices = nodes,directed = T)
net<-simplify(net, remove.loops = F)
nodes<- as_ids(V(net))


save(links,nodes,net,file = './Strict/Data/STRING_NETSCAN.RData')

load("./RowData/Databases/UniProt/uniprot_9606.tab")
