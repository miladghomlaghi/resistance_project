
load("C:/Users/kisl1/Dropbox/Resistance-CodeKarina/kinaPhosFrames_08Dic19.RData")
load("C:/Users/kisl1/Dropbox/Resistance-CodeKarina/run1to28")
library(igraph)
library(DiagrammeR)
library(DiagrammeRsvg)
library('rsvg')
library('igraph')
library(prodlim)
net_viz <- from_adj_matrix(as.matrix(get.adjacency(main_net)), mode = "directed")
png("KinasesPhos_net_test.png", units="px", width=5184, height=3546, res=600)
render_graph(net_viz,width=5184, height=3546)%>%export_svg %>% charToRaw %>% rsvg %>% png::writePNG("net.png")
dev.off()


V<-5
optionsV1 <-unlist(strsplit(final_results[V,2], ","))
range<-1:length(nodes)
all_combi<- gtools::permutations(length(nodes)-1,2, range[-V],repeats.allowed=F)
tmp<-as.numeric(optionsV1[8])
all_combi[tmp,]
combi<-c(V,all_combi[tmp,])

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

row.match(topo,minimal_match)
topo


L1_net <- delete_vertices(net, combi_names[3])
L1_short<-names(unlist(shortest_paths(L1_net,combi_names[1],combi_names[2])[[1]]))
get_sign(L1_net,combi_names[1],combi_names[2],links)
L1_short

L2_net <- delete_vertices(net, combi_names[1])
L2_short<-names(unlist(shortest_paths(L2_net,combi_names[3],combi_names[2])[[1]]))
L2_short
get_sign(L2_net,combi_names[3],combi_names[2],links)


L3_net <- delete_vertices(net, combi_names[2])
L3_short<-names(unlist(shortest_paths(L3_net,combi_names[1],combi_names[3])[[1]]))
get_sign(L3_net,combi_names[1],combi_names[3],links)
L3_short

L4_net <- delete_vertices(net, combi_names[1])
L4_short<-names(unlist(shortest_paths(L4_net,combi_names[2],combi_names[3])[[1]]))
get_sign(L4_net,combi_names[2],combi_names[3],links)
L4_short

sub_nodes<-unique(c(L1_short,L2_short,L3_short,L4_short))
sub_V5_Net<-induced_subgraph(net,sub_nodes)
plot(sub_V1_Net, vertex.size=3, edge.arrow.size=0.5,edge.curved=.3)

V(sub_V1_Net)$color = "red"



#png("V1_sub.png", units="px", width=2400, height=2400, res=300)
#dev.off()
plot_net(net,"fullNet")
plot_net(sub_V5_Net,"V5_C8Net")

[-1,]
graphAttr <- get_global_graph_attr_info(net_viz)
graphAttr <- rbind(graphAttr, c("rankdir", "LR", "graph"))

graph <- set_global_graph_attrs(graph,
                                attr = graphAttr$attr,
                                value = graphAttr$value,
                                attr_type = graphAttr$attr_type)

render_graph(graph)



plot_net<-function(net2plot,nameFile){
  net_viz <- from_adj_matrix(as.matrix(get.adjacency(net2plot)), mode = "directed")
  render_graph(net_viz,width=2400, height=2400)%>%export_svg %>% charToRaw %>% 
    rsvg %>% png::writePNG(paste(nameFile,".png",collapse = ''))
  
}






