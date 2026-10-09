
library(readxl)


data_all<-read_excel("frequency_replicates_strict.xls",col_names = T)
load("Strict/Data/topos_all.RData")
data<-data_all[data_all$Mean_Rb>0, ]

para_names<-data.frame(t(c(rep(0,9),"s","s","aa","aa","ba","ba","ca","ca","ab","ab","bb","bb","cb","cb","ac",
                           "ac","bc","bc","cc","cc", "fa","fa","fb","fb","fc","fc")))



topos_data <- data[data$`Rebound Score` > 0.3, "ID" ]

cluster_topologies <- function(topos_data, numClusters, dist_method="euclidean", hclust_method="complete"){
  
  # Ward Hierarchical Clustering
  
  d <- dist(topos_data, method=dist_method)
  
  fit <- hclust(d, method=hclust_method) 
  #fit <- hclust(d, method=method2)
  
  groups <- cutree(fit, k=numClusters) # cut tree into 5 clusters
  
  dend <- as.dendrogram(fit)
  dend <- dendextend::rotate(dend, 1:nrow(topos_data))
  colors_den = brewer.pal(8, "Dark2")
  #colGreens=colorRampPalette(c("LightGreen", "DarkMagenta")) 
  #colGreens<-filled.contour(numClusters,color.palette=colGreens)
  dend <- color_branches(dend, k=numClusters, col=colors_den)
  
  clusters = c(paste("Cluster", 1:numClusters))
  # reduce the size of the labels:
  #dend <- set(dend, "labels_cex", 0.2)
  # And plot:
  par(mar = c(1,12,1,12))
  den <- plot(dend, horiz =  TRUE, leaflab = "none")
  
  #legend("bottomleft", legend =clusters , fill = rainbow_hcl(numClusters))
  #title(main=paste(length(robustTop), " robust topologies clustered"), cex.main=1)
  
  
  
  ######## Plotting starts here#####################
  
  dend.order <- order.dendrogram(dend)
  #toposMap2<-toposMap
  #topos_data2<-topos_data*perce
  topos_data$Topology<-rownames(topos_data)
  
  levelsT = topos_data$Topology[dend.order]
  topos_data<-topos_data[dend.order,]
  toposMelt <-melt(topos_data)
  names(toposMelt)[2:3]<-c("Regulation","Value") 
  
  
  
  # Create dendrogram plot
  
  toposMelt$Topology <- factor(toposMelt$Topology, levels = levelsT,ordered = TRUE)
  
  colorLH<-c('Red4','DeepSkyblue2','White')
  
  colorM<-list(color=brewer.pal(9, "Pastel1"))
  colorM<-colorM$color[9]
  
  p<-ggplot(data = toposMelt, aes(x = Regulation, y = Topology)) +
    geom_tile(aes(fill = Value))+
    scale_fill_gradient2(name = 'Strength',low=colorLH[1],high=colorLH[2], 
                         mid=colorLH[3],breaks=c(-1,0,1),labels=c("Inhibition",
                                                                  "No connection","Activation"),
                         limits=c(-1,1)) +
    theme_linedraw(base_size = 20)+scale_y_discrete(position = "left")+
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
  
  ret_list<-list("groups"= groups, "den"= den, "cluster"= p)
  return(ret_list)
}

cluster_topologies(topos_data = topos_data, numClusters = 5)
