
source('dependencies.R')
# load all packages
#lapply(required_packages, require, character.only = TRUE)

# DATA TRANSFORMATION AND NEW VARIABLES -----------------------------------
functions_names<- list.files("./data/functions", full.names = T )
lapply(functions_names, source)

############################################################


library(igraph)
library(prodlim) #for row match
# #library(DiagrammeR)
# #library(DiagrammeRsvg)
library(foreach)
library(doParallel)
# library(readxl)
# library(tidyverse)
# library(tictoc)

source("./Functions/F07_search_exact_structures_resistance.R")
source("./functions/F19_calculate_path_value.R")
source("./functions/F04_plot_subnetwork.R")
source("./functions/F03_get_subnetwork.R")


##  importing information
# defining the nodes/proteins that should be considered in the search
start_node <- 1
end_node <- 25


data_to_match <- 'SIGNOR_NETSCAN'
topos_to_match_name <- 'topos_to_match'
num_nodes <- 3
nodesAreTargets <- c(1, 3)

# E:\postdoc\resistance\resistance_Milad\Data\matching
# load("E:/postdoc/resistance/resistance_Milad/Data/matching/SIGNOR_NETSCAN.RData")

load(paste0("./Data/matching/", data_to_match, ".RData", collapse = ''))
load(
  paste0(
    "./Data/matching/",
    num_nodes,
    "Nodes/",
    topos_to_match_name,
    ".RData",
    collapse = ''
  )
)
load(paste0("./Data/matching/", "Targets.RData", collapse = ''))
load(paste0(
  "./Data/topologies/",
  num_nodes,
  "Nodes/topos_all.RData",
  collapse = ''
))

#topos_to_match<-topos8rob_to_match
#save(topos_to_match,file='topos31_3Nodes.RData')


# finding the targetable proteins in the proposed network

targets     <- nodes[is.element(nodes, targets)]
targets_num <- match(targets, nodes)

# Modifying the variables style: changing "/" with ":"
#save(topos_all, file="topologies_4Nodes_binary")
links$ENTITYA <- gsub("/", ";", links$ENTITYA)
links$ENTITYB <- gsub("/", ";", links$ENTITYB)
links$ENTITYA <- gsub(",", ":", links$ENTITYA)
links$ENTITYB <- gsub(",", ":", links$ENTITYB)
nodes <- gsub("/", ";", nodes)
nodes <- gsub(",", ":", nodes)
targets <- gsub(",", ":", targets)
targets <- gsub("/", ";", targets)
vertex_attr(net)$name <- gsub("/", ";", vertex_attr(net)$name)
vertex_attr(net)$name <- gsub(",", ":", vertex_attr(net)$name)



###################### creating the subnetwork based on the target protein

pracma::tic()
tmp <- list()


# for (v in start_node:end_node) {
  num <- which(nodes == "AKT2")
  target_sub_graph = make_ego_graph(net,
                                    order = 2,
                                    nodes = "AKT2",
                                    mode  = c("all"))
  sub_graph_nodes <-
    igraph::as_ids(igraph::V(target_sub_graph[[1]]))
  links_sub<- as_edgelist (target_sub_graph[[1]])
  links_subnetwork <- links[row.match(as.data.frame(links_sub),links[,1:2]),]


subnodes = sub_graph_nodes
network = target_sub_graph[[1]]
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

