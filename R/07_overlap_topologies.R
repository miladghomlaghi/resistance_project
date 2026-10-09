overlap_topologies <- function(numClusters,groups,toposClus){
  
  
  # numClusters<-5
  # groups<- cluster_list$groups
  # toposClus<-topos
  # v<-"test"
  
  overlap_list<-rep(0,numClusters)
  for (i in 1:numClusters){
    #cluster <- toposClus[names(groups)[groups==i],]
    cluster <- toposClus[groups==i,]
    size_c <- nrow(cluster)
    
    overlap<-cluster[1,]
    for (j in 1:(length(cluster)-1)){
      if( sum(cluster[,j]==as.numeric(overlap[1,j]))!=nrow(cluster)){
        overlap[1,j]<-0
      }
    }
    
    #plot_binary_topology(overlap)%>%
     # export_svg %>% charToRaw %>% rsvg %>% png::writePNG(paste0("clust", i,"_",size_c,"_",v,".png", collapse =""))
    overlap_list[i]<-row.match(overlap, topos_all)
    
  
    #print(row.match(overlap, topos_all))
  }
  
  return(overlap_list)
}