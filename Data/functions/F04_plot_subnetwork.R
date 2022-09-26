plot_subnetwork<-function (subnodes, network,links, full=F){
#  source("./data/functions/F05_get_sign.R")
 # source("./data/functions/F06_get_link.R")
  
  subnodes<-unlist(subnodes)
 if(all(subnodes!='')){
   if(length(subnodes)>150){
     subnodes = subnodes[sample(1:150,20)]
     print(subnodes)
   }
  sub_net <- get_subnetwork(subnodes,network)
  sub_net <- igraph::simplify(sub_net, remove.multiple = TRUE, remove.loops = F)
 
  if(igraph::gsize(sub_net)>0){
  
  
   net <- visNetwork::toVisNetworkData(sub_net)
   
   net$edges$group<-links[prodlim::row.match(net$edges,links[,c(1,2)]) ,3]
   #for (i in 1:nrow(net$edges)){
      #net$edges$group[i] <- ifelse(net$edges$effect[i])
        #as.numeric(get_link(sub_net,net$edges[i,1],net$edges[i,2],links))
      
  # }toVisNetworkData
    
    
   net$edges$color<-ifelse(net$edges$group==1, "#62adcd","#df7f80")
    net$edges$arrows.to.type <- ifelse(net$edges$group==1,"arrow", "circle")
    net$nodes$group<-ifelse(is.element(net$nodes$id,subnodes), "N","B")
    # 
     if (full){
       color_all = "#f6f6ee"
       font_all = '18px arial black'
     } else {
       color_all = "#676b77"
       font_all = '18px arial white'
     }

    
    visNetwork(nodes = net$nodes, edges = net$edges, height = "550px",width = "600px",
               background = 'white', main = NULL)%>%
     visEdges(arrows = "to", width = 2, arrowStrikethrough = F)%>%
##7a7d8c"
      visGroups(groupname = "N", color = list(background = color_all,
                                             border = "#f6f6ee",  border = "#7a7d8c", highlight = '#e1af28'),
                shape="box",shadow = list(enabled = TRUE, size = 25), 
               font = font_all) %>%
      visGroups(groupname = "B", color = list(background = "#f6f6ee", border = "#7a7d8c",highlight = '#e1af28'),
                shape="box",shadow = list(enabled = TRUE, size = 25),
                font = '18px arial black')%>%
      visNodes(size = 80)%>%
      visPhysics(solver = "forceAtlas2Based")%>% 
     visInteraction(navigationButtons = TRUE)
  } else{
    
    sub_net <- visNetwork::toVisNetworkData(sub_net)
    visNetwork(nodes = sub_net$nodes, edges = sub_net$edges, height = "100px",width = "600px",
               background = 'white', main = list(text = "No regulations found. Select more proteins.",
                                                        style = "font-family:Tahoma;
                                                 color:#822424;font-size:20px;text-align:center;"))
   
  }}else{
    
    nodes <- data.frame(id = 1, label="First select at least two molecules in Network's proteins", color='white')
    edges <- data.frame(from = c(1,1), to = c(1,1), color='white')
    visNetwork(nodes, edges, height = "100px", width = "600px") %>% 
    visNodes(font = list(size = 22, color='#822424', family="Tahoma")) 
    
    # tmp_nodes<-igraph::V(network)
    # tmp_node<-tmp_nodes$name
    # sub_net <- igraph::delete_vertices(network, tmp_node)
    # sub_net <- visNetwork::toVisNetworkData(sub_net)
    # visNetwork(nodes = sub_net$nodes, edges = sub_net$edges, height = "10px",width = "10px",
    #            background = 'white', main = list(text = "First select at least two molecules in Network's proteins",
    #                                              style = "font-family:Tahoma; position: -50px:
    #                                              color:#822424;font-size:20px;text-align:center;"))
  }
  
}
