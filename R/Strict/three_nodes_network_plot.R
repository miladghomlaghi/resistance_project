#Figure 8, Chapter 2
#library(extrafont)
font_import()
loadfonts(device = "win")
library(DiagrammeR)
library("grDevices")
library(gr)
load("Strict/Data/topos_all.RData")


topology<- c(0, -1, 1, 0, 0, 1, 0, 1, 0)

minimal_top <- "7235  6260  6428  8039  6434  7820  7616  7841  9212  9620  9788  8009  7427  8213 12776  8258  7925  8228  13199  7460   725  6272  6863  8181  6665 15314  2841  9969  4622 10059 10581 15028 15333 14217"

minimal_top<-strsplit(minimal_top, " ")[[1]]
minimal_top<-as.numeric(minimal_top[minimal_top!=""])

top<-minimal_top[t]


top<- 8213
top_ID<-paste0("T",top, collapse = '')
topologies_3nodes<-topos_all
topology<-topologies_3nodes[top,]

A <- "<x<FONT POINT-SIZE='6'><SUB>1</SUB></FONT>>"
B<- "<x<FONT POINT-SIZE='6'><SUB>2</SUB></FONT>>"
C <- "<x<FONT POINT-SIZE='6'><SUB>3</SUB></FONT>>"

namesLinks = c("1:nw -> 1:w","2:n -> 1:e","3:n -> 1:s",
               "1:ese -> 2:wnw","2:s -> 2:ne", "3:e ->2:sw",
               "1:sw -> 3:nw", "2:w -> 3:ne", "3:w -> 3:sw")


positive <- paste(namesLinks[which(topology>0)],collapse="; ")
negative <- paste(namesLinks[which(topology<0)],collapse="; ")
null<- paste(namesLinks[which(topology==0)],collapse="; ")

#png(paste0( "Minimal",toPlot_ID[k],".png"), units="px", width=1500, height=1500, res=600)
grViz(paste0("
               digraph circles {
               
               # a 'graph' statement
               graph [overlap = true, fontsize = 11, label =",top_ID, ",fontname = 'Arial', layout = neato]
               
               node [shape = oval,fontsize = 9
               ,
               fontname = 'Arial',
               fixedsize = true, 
               style= filled, 
               color = '#676b77', fillcolor='#f6f6ee', fontcolor = '#000000',  
               height = 0.2, width = 0.4,peripheries = 0] //sets as circles
               1[label=", A," ","pos=\"0,0.8!\"]
               2[label=", B, " ","pos=\"0.75,0.4!\"]
               3[label=", C, " ","pos=\"0,0!\"]
               
               edge [color = '#62adcd',arrowsize=0.3, arrowhead=normal,  weight =2]",positive, 
             " edge [color = '#df7f80',arrowsize=0.4, arrowhead=dot,weight =2]",negative, 
             "}"))


grViz(paste0("
               digraph circles {
               
               # a 'graph' statement
               graph [overlap = true, fontsize = 2, layout = neato]
               
               node [shape = oval,fontsize = 11,
               fontname = 'Arial',
               fixedsize = true, 
               style= filled, 
               color = '#000000', fillcolor='#ffffff', fontcolor = '#000000',  
               height = 0.21, width = 0.4,peripheries = 1] //sets as circles
               1[label=", A," ","pos=\"0,0.8!\"]
               2[label=", B, " ","pos=\"0.75,0.4!\"]
               3[label=", C, " ","pos=\"0,0!\"]
               
               edge [color = '#000000',arrowsize=0.3,arrowhead=diamond,  weight =2]",positive, 
             " edge [color = '#000000',arrowsize=0.4,arrowhead=diamond, weight =2]",negative, 
             "}"))

names_vector<-c("ABL1",
                "TGFBR1",
                "MYC")
topology<-topos_all[6431,]
style_graph <- "'rounded,filled'"

  A <- names_vector[1]
  B <- names_vector[2]
  C <- names_vector[3]
  
  namesLinks = c("1:nw -> 1:w","2:n -> 1:e","3:n -> 1:s",
                 "1:ese -> 2:wnw","2:s -> 2:ne", "3:e ->2:sw",
                 "1:sw -> 3:nw", "2:w -> 3:ne", "3:w -> 3:sw")
  
  positive <- paste(namesLinks[which(topology>0)],collapse="; ")
  negative <- paste(namesLinks[which(topology<0)],collapse="; ")
  null<- paste(namesLinks[which(topology==0)],collapse="; ")
  
  #png(paste0( "Minimal",toPlot_ID[k],".png"), units="px", width=1500, height=1500, res=600)
  grViz(paste0("
               digraph circles {
               
               # a 'graph' statement
               graph [overlap = true, fontsize = 2, layout = neato]
               
               node [shape = oval,fontsize = 7,
               fontname = Helvetica,
               fixedsize = true, 
               style= filled, 
               color = '#676b77', fillcolor='#f6f6ee', fontcolor = '#676b77',  
               height = 0.2, width = 0.4,peripheries = 0] //sets as circles
               1[label=", A," ","pos=\"0,0.8!\"]
               2[label=", B, " ","pos=\"0.75,0.4!\"]
               3[label=", C, " ","pos=\"0,0!\"]
               
               edge [color = '#62adcd',arrowsize=0.3, arrowhead=normal,  weight =2]",positive, 
               " edge [color = '#df7f80',arrowsize=0.4, arrowhead=dot,weight =2]",negative, 
               "}"))

  