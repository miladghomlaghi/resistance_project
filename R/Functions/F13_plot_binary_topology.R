plot_binary_topology <-function (topology){
  
 
 if(sqrt(length(topology))==3){
   
   
 
  namesLinks = c("1:nw -> 1:w","2:n -> 1:e","3:n -> 1:s",
                 "1:se -> 2:nw","2:s -> 2:ne","3:e ->2:sw",
                 "1:sw -> 3:nw","2:w -> 3:ne","3:w -> 3:sw")
  positive <- paste(namesLinks[which(topology>0)],collapse="; ")
  negative <- paste(namesLinks[which(topology<0)],collapse="; ")
  null<- paste(namesLinks[which(topology==0)],collapse="; ")
  #png(paste0( "Minimal",toPlot_ID[k],".png"), units="px", width=1500, height=1500, res=600)
  grViz(paste0("
               digraph circles {
               
               # a 'graph' statement
               graph [overlap = true, fontsize = 2, layout = neato]
               
               node [shape = circle,fontsize = 9
               fontname = Helvetica
               fixedsize = true
               color = Gray30
               fontcolor = Gray30
               width = 0.2,peripheries = 1] //sets as circles
               1[label=x1 pos=\"0,0.8!\"]
               2[label=x2 pos=\"0.7,0.4!\"]
               3[label=x3 pos=\"0,0!\"]
               
               edge [color = SteelBlue2, arrowsize=0.3, arrowhead=normal,  weight =2]",positive, 
               " edge [color = Brown3, arrowsize=0.4, arrowhead=tee,weight =2]",negative, 
               "}"))
  
 } else if(sqrt(length(topology))== 4){
   
   
   
   namesLinks = c("1:nw -> 1:w","2:n -> 1:ne","3:n -> 1:s", "4:wnw -> 1:sse",
                  "1:e -> 2:nw","2:e -> 2:ne","3:e -> 2:sw", "4:ne -> 2:se",
                  "1:sw -> 3:nw","2:w -> 3:ne","3:w -> 3:sw", "4:w -> 3:se",
                  "1:ese -> 4:nnw","2:s -> 4:n","3:s -> 4:sw", "4:e -> 4:se")
   
   
   positive <- paste(namesLinks[which(topology>0)],collapse="; ")
   negative <- paste(namesLinks[which(topology<0)],collapse="; ")
   null<- paste(namesLinks[which(topology==0)],collapse="; ")
   #png(paste0( "Minimal",toPlot_ID[k],".png"), units="px", width=1500, height=1500, res=600)
   grViz(paste0("
               digraph circles {
               
               # a 'graph' statement
               graph [overlap = true, fontsize = 2, layout = neato]
               
               node [shape = circle
               fontsize = 9
               fontname = Helvetica
               fixedsize = true
               color = Gray30
               fontcolor = Gray30
               width = 0.2,peripheries = 1] //sets as circles
               1[label= x1 pos=\"0,1!\"]
               2[label= x2 pos=\"0.9,0.6!\"]
               3[label= x3 pos=\"0,0.4!\" ]
               4[label= x4 pos=\"0.9,0!\"]
               
               edge [color = SteelBlue2, arrowsize=0.3, arrowhead=normal,  weight =2]",positive, 
                " edge [color = Brown3, arrowsize=0.4, arrowhead=tee,weight =2]",negative, 
                "}"))
   
 } 
}