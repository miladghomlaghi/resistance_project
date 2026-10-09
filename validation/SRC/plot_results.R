
source('dependencies.R')

# DATA TRANSFORMATION AND NEW VARIABLES -----------------------------------
############################################################


library(igraph)
library(prodlim)
library(openxlsx)
library(foreach)
library(doParallel)
source("Data/functions/F03_get_subnetwork.R")
load("../../Data/initial_info.RData")


load("./validation/SRC/results_SRC.RData")


num_nodes <- 3




###################### creating the subnetwork based on the target protein

pracma::tic()


# for (v in start_node:end_node) {
# browser()

target_sub_graph = net
a<-unlist(results[c(3),])
names(a)<-NULL
unique(a)

high_robust_topos<-(results[,results[4,]>34])

# browser()
sub_graph_nodes <-high_robust_topos[1:3,1000]
high_robust_topos[1:4,1000]
links_sub<- links
links_subnetwork <- links

network = target_sub_graph

subnodes = sub_graph_nodes
sub_net = target_sub_graph
links = links_subnetwork
full = F


# browser()


subnodes<-unlist(subnodes)
if(all(subnodes!='')){
  
  sub_net <- get_subnetwork(subnodes,network)
  sub_net <- igraph::simplify(sub_net, remove.multiple = TRUE, remove.loops = F)
  
  # browser()
  if(igraph::gsize(sub_net)>0){
    
    
    net <- visNetwork::toVisNetworkData(sub_net)
    
    net$edges$group<-links[prodlim::row.match(net$edges,links[,c(1,2)]) ,3]
    
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
    
    
    visNetwork(nodes = net$nodes, edges = net$edges, height = "2000px",width = "2000px",
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
    
  }


pracma::toc()

