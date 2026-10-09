get_subnetworks_size <- function (subnetworks,network,nodes_num) {
  
  source("./R/Functions/F05_get_subnetwork.R")


  subnetworks$Size<-0
  subnetworks$NumNodes<-0
  
  for (n in 1:nrow(subnetworks)){
    sub_net <- get_subnetwork(subnetworks[n, 1:nodes_num],network, links)
    subnetworks[n,ncol(subnetworks)-1]<-gsize(sub_net)
    subnetworks[n,ncol(subnetworks)]<-length(names(V(sub_net)))
  }
  
  return(subnetworks)
  }