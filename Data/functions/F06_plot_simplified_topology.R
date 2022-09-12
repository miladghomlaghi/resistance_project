plot_simplified_topology <-function (names_vector, topology){
  
 # print(topology)
  # names_vector<-proteins
  #topology<-topo
  if(all(names_vector!='')){
  names_vector<- str_trim(stringi::stri_remove_empty(unlist(str_split(names_vector, ","))))

  
  if(sum(as.numeric(topology)!=0)>0){
    # 
    # names_vector<-subnodes
    names_vector<-gsub("-","",names_vector)
    names_vector<-gsub(" ","",names_vector)
    names_vector <- gsub(",","",names_vector)
    names_vector <- gsub("/","c",names_vector)
    
    
    if(length(names_vector)==3){
      

      namesLinks = c("1:nw -> 1:w","2:n -> 1:e","3:n -> 1:s",
                     "1:ese -> 2:wnw","2:s -> 2:ne", "3:e ->2:sw",
                     "1:sw -> 3:nw","2:w -> 3:ne","3:w -> 3:sw")
                     
                     
      
      positive <- paste(namesLinks[which(topology>0)],collapse="; ")
      negative <- paste(namesLinks[which(topology<0)],collapse="; ")
      null<- paste(namesLinks[which(topology==0)],collapse="; ")
      
      #png(paste0( "Minimal",toPlot_ID[k],".png"), units="px", width=1500, height=1500, res=600)
      grViz(paste0("
               digraph circles {
               
               # a 'graph' statement
               graph [overlap = true, fontsize = 2, layout = neato]
               
               node [shape = oval,fontsize = 6,
               fontname = Helvetica,
               fixedsize = true, 
               style= filled, 
               color = '#676b77', fillcolor='#f6f6ee', fontcolor = '#676b77',  
               height = 0.2, width = 0.4,peripheries = 0] //sets as circles
               1[label=", names_vector[1]," ","pos=\"0,0.8!\"]
               2[label=", names_vector[2], " ","pos=\"0.75,0.4!\"]
               3[label=", names_vector[3], " ","pos=\"0,0!\"]
               
               edge [color = '#62adcd',arrowsize=0.3, arrowhead=normal,  weight =2]",positive, 
                   " edge [color = '#df7f80',arrowsize=0.4, arrowhead=dot,weight =2]",negative, 
                   "}"))
      

    } else if(length(names_vector)==4){
      
      
      
      namesLinks = c("1:nw -> 1:w","2:nw -> 1:ne","3:n -> 1:s", "4:nnw -> 1:ese",
                     "1:e -> 2:w","2:e -> 2:ne","3:ene -> 2:ssw", "4:ne -> 2:se",
                     "1:sw -> 3:nw","2:wsw -> 3:nne","3:w -> 3:sw", "4:w -> 3:e",
                     "1:sse -> 4:wnw","2:s -> 4:n","3:se -> 4:sw", "4:e -> 4:se")
      #676b77
      
      positive <- paste(namesLinks[which(topology>0)],collapse="; ")
      negative <- paste(namesLinks[which(topology<0)],collapse="; ")
      null<- paste(namesLinks[which(topology==0)],collapse="; ")
      #png(paste0( "Minimal",toPlot_ID[k],".png"), units="px", width=1500, height=1500, res=600)
      grViz(paste0("
               digraph circles {
               
               # a 'graph' statement
               graph [overlap = true, fontsize = 2, layout = neato]
               
               node [shape = oval
               fontsize = 8
               fontname = Helvetica
               fixedsize = true
               color = '#f6f6ee'
               fontcolor = 'black'
               fillcolor = '#f6f6ee'
               style = filled
               height = 0.2 
               width = 0.4,peripheries = 1.1] //sets as circles
               1[label=", names_vector[1] ," pos=\"0,1!\"]
               2[label=", names_vector[2] ," pos=\"1.2,1!\"]
               3[label=", names_vector[3] ," pos=\"0,0!\" ]
               4[label=", names_vector[4] ," pos=\"1.2,0!\"]
               
               edge [color = '#62adcd', arrowsize=0.3, arrowhead=normal,  weight =2]",positive, 
                   " edge [color = '#df7f80', arrowsize=0.4, arrowhead=dot,weight =2]",negative, 
                   "}"))
      
    } else{ 
      
      main_nodes<-paste(names_vector,collapse=" ")
      
      namesLinks<-paste(permutations(length(names_vector),2,names_vector,repeats.allowed=T)[,1], "->", 
                        permutations(length(names_vector),2,names_vector,repeats.allowed=T)[,2], sep="")
      positive <- paste(namesLinks[which(topology>0)],collapse="; ")
      negative <- paste(namesLinks[which(topology<0)],collapse="; ")
      null<- paste(namesLinks[which(topology==0)],collapse="; ")
      
      
      
      grViz(paste0("
             digraph circles {
                           
                           # a 'graph' statement
                           graph [overlap = true, fontsize = 8, layout = circo]
                           
                           node [shape = oval,fontsize = 18
                           fontname = Helvetica
                           fixedsize = true,
                           height = 0.6, width = 1.1,peripheries = 1.5, color='#f6f6ee'
                           fillcolor = '#f6f6ee', style = filled, fontcolor = '#676b77', face = bold]", main_nodes," //sets as circles
                           edge [color = '#62adcd', arrowsize=1.5, arrowhead=normal,penwidth=1.5] ",positive, 
                   " edge [color = '#df7f80', arrowsize=1.5, arrowhead=dot,penwidth=1.5] ",negative, 
                   
                   
                   "}"))
                  
      
      
    }}else{
      
      plotG<- grViz("
              digraph {
              graph [label= 'No regulations found. Select more molecules.', labelloc=t,fontcolor='#822424', fontname= Tahoma,fontsize=11]
              node [shape = plaintext, fontsize = 1, fontcolor='white',fontname= Tahoma]
              A [label = '']
              B [label = '']
              C [label = '']
              D [label = '']
              E [label = '']
              edge [color = white]
              A -> B
              B -> C
              C -> D
              }"
                    )
      return(plotG)
      
    }}else{
      plotG<- grViz("
              digraph {
              graph [label= '', labelloc=t,
              fontcolor='#822424', fontname= Tahoma,fontsize=11]
              node [shape = plaintext, fontsize = 1, fontcolor='white',fontname= Tahoma]
              A [label = '']
              B [label= 'First select at least two molecules in Network`s proteins', labelloc=t,
              fontcolor='#822424', fontname= Tahoma,fontsize=7]
              C [label = '']
              edge [color = white]
              A -> B
              B -> C
              }"
      )
      return(plotG)
    }

}

