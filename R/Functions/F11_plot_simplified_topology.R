plot_simplified_topology <-function (names_vector, topology){
  
  # if (!is.vector(subnodes)){
  #   
  #   subnodes<- c(as.character(subnodes[1,1]),
  #                as.character(subnodes[1,2]),
  #                as.character(subnodes[1,3]))
  #   
  # }
  # 
  # names_vector<-subnodes
  names_vector<-gsub("-","",names_vector)
  names_vector<-gsub(" ","",names_vector)
  names_vector <- gsub(",","",names_vector)
  names_vector <- gsub("/","c",names_vector)
  
  
  
  if(length(names_vector)==3){
  
  
   A <- names_vector[1]
   B <- names_vector[2]
   C <- names_vector[3]

  namesLinks = c("1:nw -> 1:w","2:n -> 1:e","3:n -> 1:s",
                 "1:ese -> 2:wnw","2:s -> 2:ne","3:e ->2:sw",
                 "1:sw -> 3:nw","2:w -> 3:ne","3:w -> 3:sw")
  
  positive <- paste(namesLinks[which(topology>0)],collapse="; ")
  negative <- paste(namesLinks[which(topology<0)],collapse="; ")
  null<- paste(namesLinks[which(topology==0)],collapse="; ")
  #png(paste0( "Minimal",toPlot_ID[k],".png"), units="px", width=1500, height=1500, res=600)
  grViz(paste0("
               digraph circles {
               
               # a 'graph' statement
               graph [overlap = true, fontsize = 2, layout = neato]
               
               node [shape = box,fontsize = 7,
               fontname = Helvetica,
               fixedsize = true,
               style=rounded, 
               color = Gray30,  fontcolor = Gray30
               height = 0.2, width = 0.4,peripheries = 1] //sets as circles
               1[label=", A," ","pos=\"0,0.8!\"]
               2[label=", B, " ","pos=\"0.75,0.4!\"]
               3[label=", C, " ","pos=\"0,0!\"]
               
               edge [color = SteelBlue2,arrowsize=0.3, arrowhead=normal,  weight =2]",positive, 
               " edge [color = Brown3,arrowsize=0.4, arrowhead=tee,weight =2]",negative, 
               "}"))
  
  } else {
  
    A <- names_vector[1]
    B <- names_vector[2]
    C <- names_vector[3]
    D <- names_vector[4]
    
   namesLinks = c("1:nw -> 1:w","2:n -> 1:ne","3:nnw -> 1:wws", "4:nnw -> 1:se",
                  "1:e -> 2:nw","2:e -> 2:ne","3:e -> 2:sw", "4:nne -> 2:sse",
                  "1:ssw -> 3:nnw","2:w -> 3:ne","3:w -> 3:sw", "4:w -> 3:se",
                  "1:ese -> 4:nw","2:sse -> 4:nne","3:s -> 4:sw", "4:e -> 4:se")


  
  positive <- paste(namesLinks[which(topology>0)],collapse="; ")
  negative <- paste(namesLinks[which(topology<0)],collapse="; ")
  null<- paste(namesLinks[which(topology==0)],collapse="; ")
  #png(paste0( "Minimal",toPlot_ID[k],".png"), units="px", width=1500, height=1500, res=600)
  grViz(paste0("
               digraph circles {
               
               # a 'graph' statement
               graph [overlap = true, fontsize = 2, layout = neato]
               
               node [shape = box,fontsize = 7,
               fontname = Helvetica,
               fixedsize = true,
               style=rounded, 
               color = Gray30,  fontcolor = Gray30
               height = 0.2, width = 0.4,peripheries = 1] //sets as circles
               1[label=", A," ","pos=\"0,1!\"]
               2[label=", B, " ","pos=\"0.9,0.6!\"]
               3[label=", C, " ","pos=\"0,0.4!\" ]
               4[label=", D, " ","pos=\"0.9,0!\"]
               
               edge [color = SteelBlue2,arrowsize=0.3, arrowhead=normal,  weight =2]",positive, 
               " edge [color = Brown3, arrowsize=0.4, arrowhead=tee,weight =2]",negative, 
               "}"))
  }
  
}