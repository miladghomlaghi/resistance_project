plot_subnet<-function (subnodes, network,links){
  source("./R/Functions/F02_get_sign.R")
  source("./R/Functions/F03_get_link.R")

  
  sub_net <- get_subnetwork(subnodes,network, links)
  sub_net <- igraph::simplify(sub_net, remove.multiple = TRUE, remove.loops = FALSE)
  
  # if (!is.vector(subnodes)){
  #   
  #   subnodes<- c(as.character(subnodes[1,1]),
  #                as.character(subnodes[1,2]),
  #                as.character(subnodes[1,3]))
  #   }

  #sub_net<-igraph::simplify(sub_net)
  net <- toVisNetworkData(sub_net)
  net$edges$group<-0
  for (i in 1:nrow(net$edges)){
    net$edges$group[i] <- as.numeric(get_link(network,net$edges[i,1],net$edges[i,2],links))
  }
  

  subnodes<-gsub("1-phosphatidyl-1D-myo-inositol 3,4,5-trisphosphate","PI3P",subnodes)
  net$nodes$id <- gsub("1-phosphatidyl-1D-myo-inositol 3,4,5-trisphosphate","PI3P",net$nodes$id)
  net$nodes$label <- gsub("1-phosphatidyl-1D-myo-inositol 3,4,5-trisphosphate","PI3P",net$nodes$label)
  
  net$edges$from <- gsub("1-phosphatidyl-1D-myo-inositol 3,4,5-trisphosphate","PI3P",net$edges$from)
  net$edges$to <- gsub("1-phosphatidyl-1D-myo-inositol 3,4,5-trisphosphate","PI3P",net$edges$to)

  
  net$edges$color<-ifelse(net$edges$group==1, "#62adcd","#df7f80")
  net$edges$arrows.to.type <- ifelse(net$edges$group==1,"arrow", "circle")
  net$nodes$group<-ifelse(is.element(net$nodes$id,subnodes), "N","B")
  net$nodes$label<-"nnnnLLL"
  
   visNetwork(nodes = net$nodes, edges = net$edges, height = "400px",width = "400px",
              background = '#f6f6ee')%>%
     visEdges(arrows = "to", width = 5, arrowStrikethrough = F)%>%
    # visEdges(arrows = list(to = list(enabled = TRUE,
    #                                  scaleFactor = 2)))%>%
    visGroups(groupname = "N", color = list(background = "#9aa2b5", border = "#7a7d8c"),
              shape="box",shadow = list(enabled = TRUE, size = 25)) %>%
    visGroups(groupname = "B", color = list(background = "#9aa2b5", border = "#7a7d8c"),
              shape="box",shadow = list(enabled = TRUE, size = 15))%>%
    visNodes(size = 25, font = '14px arial #9aa2b5')%>%
    visPhysics(solver = "forceAtlas2Based")
   #%>% visSave(file = "./intro_network.html")
  
}