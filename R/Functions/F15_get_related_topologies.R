get_related_topologies <- function (simplified_topology, topologies){
  autoregulations<-c(1,5,9)
  regulations_values<- gtools::permutations(3,3,c(0,1,-1),repeats.allowed=T)
  
  related_topologies<-simplified_topology
  
  for (a in 1:nrow(regulations_values)){
    
    tmp_topology<-topologies[simplified_topology,]
    tmp_topology[autoregulations]<- regulations_values[a,]
    related_topologies<-c(related_topologies, prodlim::row.match(tmp_topology, as.data.frame(topologies), nomatch = 0))
  }
  
  return(unique(related_topologies))
}