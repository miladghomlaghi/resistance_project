


clustering <- function(topos_data,perce, numClusters,dist_opt, method1, method2, color_num, ret_opt,double){
  
  
  # Ward Hierarchical Clustering
  #toposMap<-toposMap*perceValue
  if (dist_opt==1){
  d <- dist(topos_data, method=method1)
 
  }else{
    distMat <- hamming.distance(data.matrix(topos_data))
    d <- as.dist(distMat)
    }# distance matrix
  fit <- hclust(d, method=method2) 
  groups <- cutree(fit, k=numClusters) # cut tree into 5 clusters
  
  dend <- as.dendrogram(fit)
  dend <- rotate(dend, 1:nrow(topos_data))
  colGreens=brewer.pal(color_num, "Dark2")
  #colGreens=colorRampPalette(c("LightGreen", "DarkMagenta")) 
  #colGreens<-filled.contour(numClusters,color.palette=colGreens)
  dend <- color_branches(dend, k=numClusters, col=colGreens)
  clusters = c(paste("Cluster", 1:numClusters))
  # reduce the size of the labels:
  #dend <- set(dend, "labels_cex", 0.2)
  # And plot:
  par(mar = c(1,12,1,12))
  plot(dend, horiz =  TRUE, leaflab = "none")
  
  #legend("bottomleft", legend =clusters , fill = rainbow_hcl(numClusters))
  #title(main=paste(length(robustTop), " robust topologies clustered"), cex.main=1)
  
  
  
  
  ######## Plotting starts here#####################
  
  dend.order <- order.dendrogram(dend)
  #toposMap2<-toposMap
  topos_data2<-topos_data*perce
  topos_data2$Topology<-rownames(topos_data)
  if(double==1){
  topos_data2[,10:18]<-(-1)*topos_data2[,10:18]}
  levelsT = topos_data2$Topology[dend.order]
  topos_data2<-topos_data2[dend.order,]
  toposMelt <-melt(topos_data2)
  names(toposMelt)[2:3]<-c("Regulation","Value") 
  
  
  
  # Create dendrogram plot
  
  toposMelt$Topology <- factor(toposMelt$Topology, levels = levelsT,ordered = TRUE)
  
  colorLH<-c('Red4','DeepSkyblue2','White')
  
  colorM<-list(color=brewer.pal(9, "Pastel1"))
  colorM<-colorM$color[9]
  
  
  
  # heatmap.plot <- ggplot(data = toposMelt, aes(x = Regulation, y = Topology)) +
  #   geom_tile(aes(fill = Value), color= "white",size=0.1) +
  #   scale_fill_gradient2(name = 'Strength',low=colorLH[1],high=colorLH[2], 
  #                        mid=colorM,breaks=c(-1,0,1),labels=c("Inhibition",
  #                                                             "No connection","Activation"),
  #                        limits=c(-1,1)) +
  #   theme_minimal(base_size = 16)+ removeGrid()+scale_y_discrete(position = "right")+
  #   theme(legend.position="left",axis.title.y=element_text(angle=90,vjust=1),axis.text.y = element_blank())
  # 
  # heatmap.plot
  
  
  p<-ggplot(data = toposMelt, aes(x = Regulation, y = Topology)) +
     geom_tile(aes(fill = Value))+
    scale_fill_gradient2(name = 'Strength',low=colorLH[1],high=colorLH[2], 
                         mid=colorLH[3],breaks=c(-1,0,1),labels=c("Inhibition",
                                                                  "No connection","Activation"),
                         limits=c(-1,1)) +
    theme_linedraw(base_size = 14)+scale_y_discrete(position = "left")+
    theme(legend.position="none",
          axis.title.y=element_text(angle=90,vjust=1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_rect(fill = "Black",
                                          colour = "Black",
                                          size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "white"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "Black"))

  
  if (ret_opt==1){
    return(groups)
  }else{
  return(p)
  }

}
plot_clust <- function(numClusters,groups,toposClust,v){
  
  
  for (i in 1:numClusters){
    cluster <- toposClus[names(groups)[groups==i],]
    size_c <- nrow(cluster)
    
    overlap<-cluster[1,]
    for (j in 1:(length(cluster)-1)){
      if( sum(cluster[,j]==as.numeric(overlap[1,j]))!=nrow(cluster)){
        overlap[1,j]<-0
      }
    }
    
    
    test1<-overlap[,1:9]
    #test1<-coreTop[4,]
    namesLinks = c("1:nw -> 1:w","2:n -> 1:e","3:n -> 1:s",
                   "1:se -> 2:nw","2:s -> 2:ne","3:e ->2:sw",
                   "1:sw -> 3:nw","2:w -> 3:ne","3:w -> 3:sw")
    positive <- paste(namesLinks[which(test1>0)],collapse="; ")
    negative <- paste(namesLinks[which(test1<0)],collapse="; ")
    null<- paste(namesLinks[which(test1==0)],collapse="; ")
    
    #png(paste("Cluster_",cl[i],".png"), units="px", width=1805, height=1200, res=300)
    grViz(paste0("
             digraph circles {
             
             # a 'graph' statement
             graph [overlap = true, fontsize = 2, layout = neato]
             
             node [shape = circle,fontsize = 9
             fontname = Helvetica
             fixedsize = true,
             width = 0.2,peripheries = 1] //sets as circles
             1[label=A pos=\"0,0.8!\"]
             2[label=B pos=\"0.7,0.4!\"]
             3[label=C pos=\"0,0!\"]
             
             edge [color = DeepSkyblue2,arrowsize=0.3, arrowhead=normal,  weight =2]",positive, 
                 " edge [color = Red4,arrowsize=0.4, arrowhead=tee,weight =2]",negative, 
                 "}"))%>%
      export_svg %>% charToRaw %>% rsvg %>% png::writePNG(paste0("clust", i,"_",size_c,"_",v,".png", collapse =""))
    
  }}
get_clustTop<- function(toposClus,groups){
  coreTop<-toposClus[1,1:9]
  for (i in 1:max(groups)){
    cluster <- toposClus[names(groups)[groups==i],]
    size_c <- nrow(cluster)
    
    overlap<-cluster[1,]
    for (j in 1:(length(cluster)-1)){
      if( sum(cluster[,j]==as.numeric(overlap[1,j]))!=nrow(cluster)){
        overlap[1,j]<-0
      }
    }
    coreTop<-rbind(coreTop,overlap[,1:9])
  }
  coreTop<-coreTop[-1,]
  return(coreTop)}



toposClus<-topos
toposClus$Topology<-rownames(toposClus)






setwd("C:/Users/kisl1/Dropbox/Resistance-CodeKarina")
#setwd("/Users/Anirka/Dropbox (Business)/Resistance-CodeKarina")
load('dataResistance.Rdata')



#namesLinks = c("1-1","2-1","3-1","1-2","2-2","3-2","1-3","2-3","3-3")
namesLinks = c("A-A","B-A","C-A","A-B","B-B","C-B","A-C","B-C","C-C")
namesLinks=as.character(namesLinks)
colnames(topos) = namesLinks




robustTop<-frec[,1]
toposMap <-topos[unlist(robustTop),] 
rownames(toposMap) <- unlist(robustTop)


robustTopPlus<-frecb[,1]
toposMapPlus <-topos[unlist(robustTopPlus),] 
rownames(toposMapPlus) <- unlist(robustTopPlus)
###########################################
##########################################
### Parameters hieracrchical clustering

#write.csv(robustTop,"robustTopAllFilter_V2_2.csv")
#write.csv(robustTopPlus,"robustTopAllFilterPlus_V2_2.csv")

groups_Rob <- clustering(toposMapPlus,perceValuePlus,5,1,"canberra","ward.D2",8,1,0)
groups_Rob
plot_clust(5,groups_Rob,toposClust,1)
coreTopIDs_Plus <-row.match(get_clustTop(toposClus,groups_Rob),topos)

top_Core_Plus <- cbind(coreTopIDs_Plus,frecuency_Rep[coreTopIDs_Plus,'Mean'])
clus_0_core <- top_Core_Plus$coreTopIDs_Plus[which(top_Core_Plus$Mean>0)]
clus_0_core_top<- topos[clus_0_core,]
row.names(clus_0_core_top) <- which(top_Core_Plus$Mean>0)
clus_0_core <- cbind(clus_0_core,which(top_Core_Plus$Mean>0))


cluster_P2 <- data.frame(cbind(robustTopPlus,groups_Rob))
names(cluster_P2)[1] <- "Topology"
tmp_list<- !is.element(cluster_P2$groups_Rob,clus_0_core[,2])
cluster_P2 <- cluster_P2[tmp_list,]
clus_P2_top<- topos[cluster_P2$Topology,]
row.names(clus_P2_top) <- cluster_P2$Topology
clus_P2_per <- perceValuePlus[tmp_list,]


groups_Rob_2 <- clustering(clus_P2_top,clus_P2_per,5,1,"canberra","ward.D2",8,2,0)
groups_Rob_2
plot_clust(6,groups_Rob_2,toposClust,11)
coreTopIDs_Plus_2 <-row.match(get_clustTop(toposClus,groups_Rob_2),topos)
top_Core_Plus_2 <- cbind(coreTopIDs_Plus_2,frecuency_Rep[coreTopIDs_Plus_2,'Mean'])
clus_02_core<- top_Core_Plus_2$coreTopIDs_Plus_2[which(top_Core_Plus_2$Mean>0)]
clus_02_core <- data.frame(cbind(clus_02_core,which(top_Core_Plus_2$Mean>0),2))
names(clus_02_core)<-c('Topology','Cluster', 'NumCluster')




ward.D2

perceData<- cbind(perceValue,perceValue)
toposData<- toposMap2


perceData<- (perceValue-perceValue)+1
toposData<- toposMap

png("TreeClus1.png", units="px", width=800, height=1200, res=300)
groups <- clustering(toposData,perceData,100,1,"canberra","ward.D2",8,2,0)
png("HeatClus1.png", units="px", width=1200, height=1200, res=300)
groups
dev.off()
#plot_clust(200,groups,toposClust,1)
coreTopIDs <-row.match(get_clustTop(toposClus,groups),topos)
top_Core <- cbind(coreTopIDs,frecuency_Rep[coreTopIDs,'Mean'])
clus_1_core <- top_Core$coreTopIDs[which(top_Core$Mean>0)]
clus_1_core <- data.frame(cbind(clus_1_core,which(top_Core$Mean>0),1))
names(clus_1_core)<-c('Topology','Cluster', 'NumCluster')
#cluster<- topos[unlist(robustTop),] 
#cluster$Topology<- 1:nrow(robustTop)


#### Best clusters with ward.D2 and canberra

cluster_2 <- cbind(robustTop,groups)
tmp_list<- !is.element(cluster_2$groups,clus_1_core[,2])
cluster_2 <- cluster_2[tmp_list,]
clus_2_top <- topos[cluster_2$Topology,]
row.names(clus_2_top) <- cluster_2$Topology
clus_2_per <- perceData[tmp_list,]


group_2 <- clustering(clus_2_top,clus_2_per,100,1,"canberra","ward.D2",8,2,0)
png("HeatClus2.png", units="px", width=1200, height=1200, res=300)
group_2
dev.off()
#plot_clust(100,group_2,toposClust,5)
coreTopIDs_2 <-row.match(get_clustTop(toposClus,group_2),topos)
top_Core_2 <- cbind(coreTopIDs_2,frecuency_Rep[coreTopIDs_2,'Mean'])
clus_2_core <- top_Core_2$coreTopIDs_2[which(top_Core_2$Mean>0)]
clus_2_core <- data.frame(cbind(clus_2_core,which(top_Core_2$Mean>0),2))
names(clus_2_core)<-c('Topology','Cluster', 'NumCluster')


cluster_3 <- data.frame(cbind(cluster_2$Topology,group_2))
names(cluster_3)[1] <- "Topology"
tmp_list<- !is.element(cluster_3$group_2,clus_2_core[,2])
cluster_3 <- cluster_3[tmp_list,]
clus_3_top<- topos[cluster_3$Topology,]
row.names(clus_3_top) <- cluster_3$Topology
clus_3_per <- clus_2_per[tmp_list,]


group_3 <- clustering(clus_3_top,clus_3_per,100,1,"canberra","ward.D2",8,2,0)
png("HeatClus3.png", units="px", width=1200, height=1200, res=300)
group_3
dev.off()
#plot_clust(48,group_3,toposClust,6)
coreTopIDs_3 <-row.match(get_clustTop(toposClus,group_3),topos)
match(coreTopIDs_3,top_Core_2$coreTopIDs_2)
top_Core_3 <- cbind(coreTopIDs_3,frecuency_Rep[coreTopIDs_3,'Mean'])
clus_3_core <- top_Core_3$coreTopIDs_3[which(top_Core_3$Mean>0)]
clus_3_core <- data.frame(cbind(clus_3_core,which(top_Core_3$Mean>0),3))
names(clus_3_core)<-c('Topology','Cluster', 'NumCluster')



cluster_4 <- data.frame(cbind(cluster_3$Topology,group_3))
names(cluster_4)[1] <- "Topology"
tmp_list<- !is.element(cluster_4$group_3,clus_3_core[,2])
cluster_4 <- cluster_4[tmp_list,]
clus_4_top<- topos[cluster_4$Topology,]
row.names(clus_4_top) <- cluster_4$Topology
clus_4_per <- clus_3_per[tmp_list,]


group_4 <- clustering(clus_4_top,clus_4_per,100,1,"canberra","ward.D2",8,2,0)
png("HeatClus4.png", units="px", width=1200, height=1200, res=300)
group_4
dev.off()
#plot_clust(48,group_3,toposClust,6)
coreTopIDs_4 <-row.match(get_clustTop(toposClus,group_4),topos)
match(coreTopIDs_4,top_Core_3$coreTopIDs_3)
top_Core_4 <- cbind(coreTopIDs_4,frecuency_Rep[coreTopIDs_4,'Mean'])
clus_4_core <- top_Core_4$coreTopIDs_4[which(top_Core_4$Mean>0)]
clus_4_core <- data.frame(cbind(clus_4_core,which(top_Core_4$Mean>0),4))
names(clus_4_core)<-c('Topology','Cluster', 'NumCluster')


cluster_5 <- data.frame(cbind(cluster_4$Topology,group_4))
names(cluster_5)[1] <- "Topology"
tmp_list<- !is.element(cluster_5$group_4,clus_4_core[,2])
cluster_5 <- cluster_5[tmp_list,]
clus_5_top<- topos[cluster_5$Topology,]
row.names(clus_5_top) <- cluster_5$Topology
clus_5_per <- clus_4_per[tmp_list,]

group_5 <- clustering(clus_5_top,clus_5_per,100,1,"canberra","ward.D2",8,2,0)
png("HeatClus5.png", units="px", width=1200, height=1200, res=300)
group_5
dev.off()

#plot_clust(48,group_5,toposClust,8)
coreTopIDs_5 <-row.match(get_clustTop(toposClus,group_5),topos)
top_Core_5 <- cbind(coreTopIDs_5,frecuency_Rep[coreTopIDs_5,'Mean'])
clus_5_core <- top_Core_5$coreTopIDs_5[which(top_Core_5$Mean>0)]
clus_5_core <- data.frame(cbind(clus_5_core,which(top_Core_5$Mean>0),5))
names(clus_5_core)<-c('Topology','Cluster', 'NumCluster')

# 
 cluster_6 <- data.frame(cbind(cluster_5$Topology,group_5))
 names(cluster_6)[1] <- "Topology"
 tmp_list<- !is.element(cluster_6$group_5,clus_5_core[,2])
 cluster_6 <- cluster_6[tmp_list,]
 clus_6_top<- topos[cluster_6$Topology,]
 row.names(clus_6_top) <- cluster_6$Topology
 clus_6_per <- clus_5_per[tmp_list,]
# 



group_6 <- clustering(clus_6_top,clus_6_per,75,1,"canberra","ward.D2",8,2,0)
png("HeatClus6.png", units="px", width=1200, height=1200, res=300)
group_6
dev.off()


# plot_clust(48,group_6,toposClust,9)
 coreTopIDs_6 <-row.match(get_clustTop(toposClus,group_6),topos)
 top_Core_6 <- cbind(coreTopIDs_6,frecuency_Rep[coreTopIDs_6,'Mean'])
 clus_6_core <- top_Core_6$coreTopIDs_6[which(top_Core_6$Mean>0)]
 clus_6_core <- data.frame(cbind(clus_6_core,which(top_Core_6$Mean>0),6))
 names(clus_6_core)<-c('Topology','Cluster', 'NumCluster')
# 
cluster_7 <- data.frame(cbind(cluster_6$Topology,group_6))
names(cluster_7)[1] <- "Topology"
tmp_list<- !is.element(cluster_7$group_6,clus_6_core[,2])
cluster_7 <- cluster_7[tmp_list,]
clus_7_top<- topos[cluster_7$Topology,]
row.names(clus_7_top) <- cluster_7$Topology
clus_7_per <- clus_6_per[tmp_list,]

# 
clus_7_core <-data.frame(cbind(cluster_7$Topology,1:nrow(cluster_7),7))
names(clus_7_core)<-c('Topology','Cluster', 'NumCluster')
#

colorLH<-c('Red4','DeepSkyblue2','White')
clust_7_topH<-clus_7_top
clust_7_topH$Topology<-rownames(clust_7_topH)
toposMelt7 <-melt(clust_7_topH)
names(toposMelt7)[2:3]<-c("Regulation","Value") 
png("HeatClus7.png", units="px", width=1200, height=1200, res=300)
ggplot(data = toposMelt7, aes(x = Regulation, y = Topology)) +
  geom_tile(aes(fill = Value)) +
  scale_fill_gradient2(name = 'Strength',low=colorLH[1],high=colorLH[2], 
                       mid=colorLH[3],breaks=c(-1,0,1),labels=c("Inhibition",
                                                                "No connection","Activation"),
                       limits=c(-1,1)) +
  theme_linedraw(base_size = 14)+scale_y_discrete(position = "left")+
  theme(legend.position="none",
        axis.title.y=element_text(angle=90,vjust=1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = "White",
                                        colour = "Black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "Black"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "Black"))
dev.off()


install.packages('dbscan')
library("dbscan")

cl <- hdbscan(toposMap, minPts = 2)

data("moons")
