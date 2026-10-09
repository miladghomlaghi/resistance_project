setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load('R_Data/kinaPhosFrames_10Dic19.RData')
load("frecuencyMatchDesk")
load("matching_scoresNames")
gene_result <- read.delim("C:/Users/kisl1/Downloads/PhD/Resistance-CodeKarina/R_codes/R_Data/gene_result.txt")
gene_result_sup <- read.delim("C:/Users/kisl1/Downloads/PhD/Resistance-CodeKarina/R_codes/R_Data/gene_result_sup.txt")

matching_scores<-cbind(frecuency_match[,1], matching_scores)
matching_scores<-matching_scores[-c(39,47,48),]
final_Frequ_match <- cbind(frecuency_match[,1],rowSums(frecuency_match[,2:16]))



toPlot<-topos_to_match
toPlot_ID<-row.match(topos_to_match,topos)
frecuency_match<- data.frame(toPlot_ID)

matching_Datafiles <- list.files(path="R_Data/TopoMatchData")
for (f in 1:length(matching_Datafiles)){
  load(paste0('R_Data/TopoMatchData/',matching_Datafiles[f],collapse =''))
  frecuency_match[,1+f]<-0
  
  for (V in 1:nrow(final_results)){
    optionsV1 <-unlist(strsplit(final_results[V,2], ","))
    if (optionsV1 != "none"){
      range<-1:length(nodes)
      vv<-as.numeric(final_results[V,1])
      all_combi<- gtools::permutations(length(nodes)-1,2, range[-vv],repeats.allowed=F) 
      for (C in 1:length(optionsV1)){
        tmp<-as.numeric(optionsV1[C])
        combi<-c(vv,all_combi[tmp,])
        
        combi_names<-nodes[combi]
        topo <-rep(0,9)
        tmp_net <- delete_vertices(net, combi_names[3])
        topo[2] <- get_link(tmp_net,combi_names[2],combi_names[1],links)
        tmp_net <- delete_vertices(net, combi_names[2])
        topo[3] <- get_link(tmp_net,combi_names[3],combi_names[1],links)
        tmp_net <- delete_vertices(net, combi_names[3])
        topo[4] <- get_link(tmp_net,combi_names[1],combi_names[2],links)
        tmp_net <- delete_vertices(net, combi_names[1])
        topo[6] <- get_link(tmp_net,combi_names[3],combi_names[2],links)
        tmp_net <- delete_vertices(net, combi_names[2])
        topo[7] <- get_link(tmp_net,combi_names[1],combi_names[3],links)
        tmp_net <- delete_vertices(net, combi_names[1])
        topo[8] <- get_link(tmp_net,combi_names[2],combi_names[3],links)
        
        
        frecuency_match[row.match(topo,topos_to_match),1+f]<-frecuency_match[row.match(topo,topos_to_match),1+f]+1
      }}}
}

save(frecuency_match,file = "frecuencyMatchDesk")



combi_names<- c("IRS1",
                "MTOR",
                "AKT1")

nodes_combi_split_filter2<-nodes_combi_perGene
nodes_combi_split_filter2<-nodes_combi_split
nodes_combi_split_filter2<-cbind(nodes_combi_split_filter2,numEdges=0) #add extra column
nodes_combi_split_filter2<-cbind(nodes_combi_split_filter2,numNodes=0) #add extra column

 netPlot<-netStringSignorKinase
 linksPlot<-linksStringSignorKinase
for (h in 1:nrow(nodes_combi_split_filter2)){
  

# a -----------------------------------------------------------------------

  
  h<-small_matches[h]
  combi_names<- nodes_combi_split_filter2[h,]
  
 
  L1_net <- delete_vertices(netPlot, combi_names[3])
  L1_short<-names(unlist(shortest_paths(L1_net,combi_names[1],combi_names[2])[[1]]))
  L1_shortB<-names(unlist(shortest_paths(L1_net,combi_names[2],combi_names[1])[[1]]))
  L1_short
  L1_shortB
  get_sign(L1_net,combi_names[1],combi_names[2],linksPlot)
  get_sign(L1_net,combi_names[2],combi_names[1],linksPlot)
  
  L2_net <- delete_vertices(netPlot, combi_names[2])
  L2_short<-names(unlist(shortest_paths(L2_net,combi_names[1],combi_names[3])[[1]]))
  L2_shortB<-names(unlist(shortest_paths(L2_net,combi_names[3],combi_names[1])[[1]]))
  L2_short
  L2_shortB
  #get_sign(L2_net,combi_names[1],combi_names[3],linksPlot)
  
  
  L3_net <- delete_vertices(netPlot, combi_names[1])
  L3_short<-names(unlist(shortest_paths(L3_net,combi_names[2],combi_names[3])[[1]]))
  L3_shortB<-names(unlist(shortest_paths(L3_net,combi_names[3],combi_names[2])[[1]]))
  get_sign(L3_net,combi_names[2],combi_names[3],linksPlot)
  get_sign(L3_net,combi_names[3],combi_names[2],linksPlot)
  L3_short
  L3_shortB
  
  sub_nodes<-unique(c(L1_short,L1_shortB,L2_short,L2_shortB,L3_short,L3_shortB))
  sub_V1_Net<-induced_subgraph(netPlot,sub_nodes)
  
  net_Test <- toVisNetworkData(sub_V1_Net)
  net_Test$edges$group<-0
  for (i in 1:nrow(net_Test$edges)){
    net_Test$edges$group[i] <- as.numeric(get_link(netPlot,net_Test$edges[i,1],net_Test$edges[i,2],linksPlot))
  }
  
  net_Test$edges$color<-ifelse(net_Test$edges$group==1, "#00b2ee","#8b0000")
  net_Test$nodes$group<-ifelse(is.element(net_Test$nodes$id,combi_names), "N","B")
  
  
  
  visNetwork(nodes = net_Test$nodes, edges = net_Test$edges, height = "600px",width = "600px")%>%
    visEdges(arrows = list(to = list(enabled = TRUE, 
                                     scaleFactor = 1.5, type = "arrow")))%>% 
    visGroups(groupname = "N", color = list(background = "lightgrey", border = "black"), shape="box",shadow = list(enabled = TRUE, size = 25)) %>%
    visGroups(groupname = "B", color = list(background = "white", border = "grey"),shape="box",shadow = list(enabled = TRUE, size = 15))%>%
    visNodes(size = 25,font = list(size = 16))%>%
    visPhysics(solver = "forceAtlas2Based")
  
  

# b -----------------------------------------------------------------------

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #sub_V1_Net<-simplify(sub_V1_Net)
                       #, remove.multiple = TRUE, remove.loops = TRUE,
                    #edge.attr.comb = igraph_opt("edge.attr.comb"))
  
  
  temTop<-topos[1,]
  temTop[1]<-0
  temTop[2]<-get_link(L1_net,combi_names[2],combi_names[1],linksPlot)
  temTop[3]<-get_link(L2_net,combi_names[3],combi_names[1],linksPlot)
  temTop[4]<-get_link(L1_net,combi_names[1],combi_names[2],linksPlot)
  temTop[5]<-0
  temTop[6]<-get_link(L3_net,combi_names[3],combi_names[2],linksPlot)
  temTop[7]<-get_link(L2_net,combi_names[1],combi_names[3],linksPlot)
  temTop[8]<-get_link(L3_net,combi_names[2],combi_names[3],linksPlot)
  temTop[9]<-0
  
  
  nodes_combi_split_filter2[h,4]<-gsize(sub_V1_Net)
  nodes_combi_split_filter2[h,5]<-length(names(V(sub_V1_Net)))
 # plot(sub_V1_Net)
}

save(nodes_combi_split_filter2,file = "combinationsStrSigKin8213")

nodes_combi_split_filter2DF<-data.frame(nodes_combi_split_filter2)
nodes_combi_split_filter2DF$numNodes<-as.numeric(as.character(nodes_combi_split_filter2DF$numNodes))
nodes_combi_split_filter2DF$numEdges<-as.numeric(as.character(nodes_combi_split_filter2DF$numEdges))
min(nodes_combi_split_filter2DF$numNodes)
small_matches<-which(nodes_combi_split_filter2DF$numNodes<6)

#V(sub_V1_Net)$color<-ifelse(is.element(names(V(sub_V1_Net)),combi_names), "black","grey")
#E(sub_V1_Net)$color<-ifelse(net_Test$edges$group==1, "#00b2ee","#8b0000")

# plot.igraph(sub_V1_Net,edge.curved=.2, layout=layout_nicely,vertex.shape="none",
#             vertex.label.color=V(sub_V1_Net)$color, vertex.label.font=1,
#             vertex.label.cex=1.3, edge.width=2,edge.arrow.size=.6,
#             edge.arrow.width=1.2, edge.label.family="Heveltica",edge.arrow.size=1.5)
# 
# net_viz <- from_adj_matrix(as.matrix(get.adjacency(sub_V1_Net)), mode = "directed")

net_Test <- toVisNetworkData(sub_V1_Net)
net_Test$edges$group<-0
for (i in 1:nrow(net_Test$edges)){
  net_Test$edges$group[i] <- as.numeric(get_link(netPlot,net_Test$edges[i,1],net_Test$edges[i,2],linksPlot))
}

net_Test$edges$color<-ifelse(net_Test$edges$group==1, "#00b2ee","#8b0000")
net_Test$nodes$group<-ifelse(is.element(net_Test$nodes$id,combi_names), "N","B")



visNetwork(nodes = net_Test$nodes, edges = net_Test$edges, height = "600px",width = "600px")%>%
  visEdges(arrows = list(to = list(enabled = TRUE, 
                                   scaleFactor = 1.5, type = "arrow")))%>% 
  visGroups(groupname = "N", color = list(background = "lightgrey", border = "black"), shape="box",shadow = list(enabled = TRUE, size = 25)) %>%
  visGroups(groupname = "B", color = list(background = "white", border = "grey"),shape="box",shadow = list(enabled = TRUE, size = 15))%>%
  visNodes(size = 25,font = list(size = 16))%>%
  visPhysics(solver = "forceAtlas2Based")



visNetwork(nodes = net_Test$nodes, edges = net_Test$edges, height = "600px",width = "600px")%>%
  visEdges(arrows = list(to = list(enabled = TRUE, 
                                   scaleFactor = 1.5, type = "arrow")))%>% visGroups(groupname = "N", color = list(background = "lightgrey", border = "black", borderWidth = 7), shape="box") %>%
  visGroups(groupname = "B", color = list(background = "white", border = "grey"),shape="box")%>%
  visNodes(shadow = list(enabled = TRUE, size = 25), size = 25,font = list(size = 16))



visGroups(groupname = "N", color = "red", shape="text") %>%
  visGroups(groupname = "B", color = "grey", shape="text")%>% 
  visNodes(shadow = list(enabled = TRUE, size = 20), size = 25)



visGroups(groupname = "N", color = list(background = "white", border = "black")) %>%
  visGroups(groupname = "B", color = list(background = "white", border = "grey",borderWidth = 1),shape="text")%>% visNodes(shadow = list(enabled = TRUE, size = 20), size = 25)



visEdges(arrows ="to")
%>%                            # arrow "to" for all edges
  visGroups(groupname = "A", color = "darkblue") %>%    # darkblue for group "A"
  visGroups(groupname = "B", color = "red")       
visIgraph(sub_V1_Net)
render_graph(net_viz)

# png("KinasesPhos_net.png", units="px", width=2400, height=2400, res=300)
# render_graph(net_viz,width=2400, height=2400)%>%export_svg %>% charToRaw %>% rsvg %>% png::writePNG("net.png", collapse ="")
# dev.off()

final_Frequ_match<-frecuency_match
types <- c("PFBL",
           "NFBL",
           "DNFBL",
           "SNFB",
           "IFFL",
           "CFFL",
           "NFBL+NFBL",
           "NFBL+IFFL",
           "NFBL+DNFBL",
           "IFFL+PFBL",
           "NFBL+PFBL"
)
classtest<-c(rep("a",10),rep("b",10),rep("c",11))

# Create dataset
data <- data.frame(
  individual=paste( "Mister ", seq(1,60), sep=""),
  group=c( rep('A', 10), rep('B', 30), rep('C', 14), rep('D', 6)) ,
  value=sample( seq(10,100), 60, replace=T)
)

class<-c(types[2],types[2],types[5],types[2])

final_Frequ_match0<-final_Frequ_match[final_Frequ_match[,2]>0,]
final_Frequ_match0<-final_Frequ_match0[-15,]
data <- data.frame(
  individual=as.factor(final_Frequ_match0[,1]),
  group=class,
  value=(final_Frequ_match0[,2]*100)/max(final_Frequ_match0[,2])
)
data = data %>% arrange(group, value)


# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 2
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))


empty_bar <- 2
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))



# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)
# prepare a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Make the plot
p <- ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  #annotate("text", x = rep(max(data$value),4), y = c(20, 40, 60, 80), label = c("20", "40", "60", "80") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.6) +
  ylim(-100,120) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=-26, label=individual, hjust=hjust, color=group),  fontface="bold",alpha=0.9, size=1.4, angle= label_data$angle, inherit.aes = FALSE ) +
  geom_text(data=data, aes(x=id, y=value+5, label=value*max(final_Frequ_match0[,2])/100, hjust=label_data$hjust),color="black", fontface="bold",alpha=0.9, size=1.4, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start-0.5, y = -29, xend = end+0.5, yend = -29, color=group),  alpha=0.6, size=0.5 , inherit.aes = FALSE )  
+
  geom_text(data=base_data, aes(x = title, y =-30,label=group),hjust=c(1,1,1,0,0,0), colour = "black", alpha=0.6, size=2,angle= c(20,45,-60,-60,-20,0), fontface="bold", inherit.aes = FALSE)

png("ciruclarFreqV2.png", units="px", width=1600, height=2400, res=600)
p
dev.off()



### All the combinations for a specific Top

load("R_Data/matching_scoresNamesStrSigKinase")
matching_scores<-frecuency_match
top<- 8213
nodes_combi<- matching_scores[matching_scores[,1]==top,3]

nodes_combi_split <- c("Node 1","Node 2", "Node 3")
  tmp <- nodes_combi
  tmp_split <- unlist(strsplit(tmp,"/"))
  for (s in 2: length(tmp_split)){
    tmp2 <- tmp_split[s]
    tmp2_split<-unlist(strsplit(tmp2,","))
    nodes_combi_split <- rbind(nodes_combi_split,tmp2_split)
}

colnames(nodes_combi_split)<-c("Node 1","Node 2", "Node 3")
nodes_combi_split<-nodes_combi_split[-1,]
row.names(nodes_combi_split)<-1:nrow(nodes_combi_split)


nodes_combi_perGene<-nodes_combi_split[nodes_combi_split[,1]=="MAP4K5",]
rownames(nodes_combi_perGene)<-1:nrow(nodes_combi_perGene)
nodes_combi_perGene<-nodes_combi_perGene[1:60,]

nodes_combi_split_filter<-nodes_combi_split[1,]
for (g in 1:nrow(nodes_combi_split)){
  
  node_name <- nodes_combi_split[g,3]
  if(length(which(gene_result$Symbol==node_name))){
    nodes_combi_split_filter<-rbind(nodes_combi_split_filter,nodes_combi_split[g,])
  }
}
nodes_combi_split_filter<-nodes_combi_split_filter[-1,]


nodes_combi_split_filter2<-nodes_combi_split[1,]
for (g in 1:nrow(nodes_combi_split_filter)){
  
  node_name <- nodes_combi_split_filter[g,1]
  if(length(which(gene_result$Symbol==node_name))){
    nodes_combi_split_filter2<-rbind(nodes_combi_split_filter2,nodes_combi_split_filter[g,])
  }
}
nodes_combi_split_filter2<-nodes_combi_split_filter2[-1,]
row.names(nodes_combi_split_filter2)<-1:nrow(nodes_combi_split_filter2)
