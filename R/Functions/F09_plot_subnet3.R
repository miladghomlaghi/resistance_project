plot_subnet3<-function (subnodes, net,links, save=F){
  source("./R/Functions/F02_get_sign.R")
  source("./R/Functions/F03_get_link.R")
  

  sub_net <- get_subnetwork(subnodes,net)
  sub_net <- igraph::simplify(sub_net, remove.multiple = TRUE, remove.loops = FALSE)
  
  # if (!is.vector(subnodes)){
  #   
  #   subnodes<- c(as.character(subnodes[1,1]),
  #                as.character(subnodes[1,2]),
  #                as.character(subnodes[1,3]))
  #   
  # }
  
  netB <- toVisNetworkData(sub_net)
  netB$edges$group<-0
  for (i in 1:nrow(netB$edges)){
    netB$edges$group[i] <- as.numeric(get_link(net,netB$edges[i,1],netB$edges[i,2],links))
  } 

#a<- as_data_frame(netB, what = c("edges", "nodes", "both"))
a<- netB$edges
a
names_vector<- unique(c(a[,1],a[,2]))

subnodes<-gsub("1-phosphatidyl-1D-myo-inositol 3,4,5-trisphosphate","PI3P",subnodes)
names_vector <- gsub("1-phosphatidyl-1D-myo-inositol 3,4,5-trisphosphate","PI3P",names_vector)
a$from <- gsub("1-phosphatidyl-1D-myo-inositol 3,4,5-trisphosphate","PI3P",a$from)
a$to <- gsub("1-phosphatidyl-1D-myo-inositol 3,4,5-trisphosphate","PI3P",a$to)

subnodes<-gsub(" ","",subnodes)
names_vector <- gsub(" ","",names_vector)
a$from <- gsub(" ","",a$from)
a$to <- gsub(" ","",a$to)

subnodes<-gsub("-","",subnodes)
names_vector <- gsub("-","",names_vector)
a$from <- gsub("-","",a$from)
a$to <- gsub("-","",a$to)

subnodes<-gsub(",","",subnodes)
names_vector <- gsub(",","",names_vector)
a$from <- gsub(",","",a$from)
a$to <- gsub(",","",a$to)

subnodes<-gsub("/","c",subnodes)
names_vector <- gsub("/","c",names_vector)
a$from <- gsub("/","c",a$from)
a$to <- gsub("/","c",a$to)



# A <- subnodes[1]
# B <- subnodes[2]
# C <- subnodes[3]
main_nodes<-paste(subnodes,collapse=" ")
no_main <- paste(names_vector[!is.element(names_vector,subnodes)],collapse=" ")
namesLinks = paste0(a[,1],"->",a[,2])


positive <- paste(namesLinks[which(netB$edges$group>0)],collapse="; ")
negative <- paste(namesLinks[which(netB$edges$group<0)],collapse="; ")
#null<- paste(namesLinks[which(test1==0)],collapse="; ")
#png(paste0( "Minimal",toPlot_ID[k],".png"), units="px", width=1500, height=1500, res=600)
#png(paste0(paste0(subnodes,collapse = "-"),".png"), units="px", width=800, height=1200, res=300)

if(save){
  grViz(paste0("
             digraph circles {
             
             # a 'graph' statement
             graph [overlap = true, fontsize = 8, layout = dot]
             
             node [shape = oval,fontsize = 18
             fontname = Helvetica
             fixedsize = true,
             height = 0.6, width = 1.1,peripheries = 1.5] //sets as circles
             node [fillcolor = #9aa2b5, style = filled] ",main_nodes,
               " node [fillcolor = #f6f6ee]", no_main,
               " edge [color = #62adcd,arrowsize=1, arrowhead=normal,penwidth=1.5] ",positive, 
               " edge [color = #df7f80, arrowsize=1.5, arrowhead=tee,penwidth=1.5] ",negative, 
               "}"))%>%export_svg %>% charToRaw %>% rsvg %>% png::writePNG(
                 paste0("./images/",paste0(subnodes,collapse = "-"),".png"), dpi =300)
}else{
  
  plotG<-  grViz(paste0("
             digraph circles {
                           
                           # a 'graph' statement
                           graph [overlap = true, fontsize = 8, layout = dot]
                           
                           node [shape = oval,fontsize = 18
                           fontname = Helvetica
                           fixedsize = true,
                           height = 0.6, width = 1.1,peripheries = 1.5] //sets as circles
                           node [fillcolor = '#9aa2b5'', style = filled] ",main_nodes,
                        " node [fillcolor = f6f6ee]", no_main,
                        " edge [color = 62adcd,arrowsize=1, arrowhead=normal,penwidth=1.5] ",positive, 
                        " edge [color = df7f80,arrowsize=1.5, arrowhead=tee,penwidth=1.5] ",negative, 
                        "}"))
  plotG
  return(plotG)
    }                                                                    
}

