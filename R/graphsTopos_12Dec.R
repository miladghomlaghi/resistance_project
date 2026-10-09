



# random graph
net = rgraph(3, mode = "digraph", tprob = 0.5)
net = network(net, directed = FALSE)

# vertex names
network.vertex.names(net) = letters[1:3]

ggnet2(net)
test1<- toposMap[,-10][1,]

test1<-overlap[,1:9]
topNum<-8213
test1 <- topos[topNum,]


test1<-temTop
#topos_ID <-results_UU$Topology
#toPlot <- topos[topos_ID,]
topos_ID <- frecb$Topology
toPlot<-robustTop
toPlot<-topos_to_match
toPlot<-minimal_match
toPlot_ID<-row.match(minimal_match,topos)
for (k in 1:nrow(toPlot)){
  
  test1<-toPlot[k,]
  
namesLinks = c("1:nw -> 1:w","2:n -> 1:e","3:n -> 1:s",
               "1:se -> 2:nw","2:s -> 2:ne","3:e ->2:sw",
               "1:sw -> 3:nw","2:w -> 3:ne","3:w -> 3:sw")
positive <- paste(namesLinks[which(test1>0)],collapse="; ")
negative <- paste(namesLinks[which(test1<0)],collapse="; ")
null<- paste(namesLinks[which(test1==0)],collapse="; ")
#png(paste0( "Minimal",toPlot_ID[k],".png"), units="px", width=1500, height=1500, res=600)
grViz(paste0("
digraph circles {

  # a 'graph' statement
  graph [overlap = true, fontsize = 2, layout = neato]

  node [shape = circle,fontsize = 9
  fontname = Helvetica
        fixedsize = true,
        width = 0.2,peripheries = 1] //sets as circles
  1[label=A pos=\"0,0.8!\"]
  2[label=B pos=\"0.7,0.4!\"]
  3[label=C pos=\"0,0!\"]

  edge [color = DeepSkyblue2,arrowsize=0.3, arrowhead=normal,  weight =2]",positive, 
  " edge [color = Red4,arrowsize=0.4, arrowhead=tee,weight =2]",negative, 
"}"))


%>%export_svg %>% charToRaw %>% rsvg %>% png::writePNG(paste0( "Minimal",topNum,".png"))


dev.off()
  
}
  
  
  
  export_svg %>% charToRaw %>% rsvg %>% png::writePNG(paste0( "Minimal",toPlot_ID[k],".png", collapse =""),units="px", width=1500, height=1500, res=600)
}




grViz(paste0("
digraph circles {

  # a 'graph' statement
  graph [overlap = true, fontsize = 2, layout = neato]

  node [shape = circle,fontsize = 9
  fontname = Helvetica
        fixedsize = true,
        width = 0.2,peripheries = 1] //sets as circles
  1[label=A pos=\"0,0.8!\"]
  2[label=B pos=\"0.7,0.4!\"]
  3[label=C pos=\"0,0!\"]

  edge [color = black,arrowsize=0.3, arrowhead=vee,  weight =2]",positive, 
             " edge [color = black,arrowsize=0.4, arrowhead=vee,weight =2]",negative, 
             "}"))













colors <- c("CornflowerBlue","BlueViolet")

ggplot(over, aes(x=Mean_OV, y=Mean_Ratio)) + geom_point()+theme_minimal()

png("over_1.png", units="px", width=1805, height=1200, res=300)
#ggplot(over, aes(x=Mean_OV, y=Mean_BB)) + geom_point(color=ifelse(over$Mean_BB < 25,"CornflowerBlue","BlueViolet"))+theme_minimal()+xlim(c(0,120))+ylim(c(0,120))
ggplot(over, aes(x=Mean_BB, y=Mean_OV)) + 
  geom_point(color=ifelse(over$Mean < 0.01,"DarkMagenta","DarkOrange"), size=0.9)+
  theme_minimal()+xlim(c(0,120))+ylim(c(0,120))+
   xlab("Bounce-Back sets") + ylab("Overshoot sets")+ 
  geom_abline(intercept=0, slope = 1,linetype='dashed', color='black', size=0.8)

dev.off()

"DarkMagenta","DarkOrange"

over$Class <-"Normal"
over[over$Mean>0.01,'Class']<-"Robust"
over$Mean <-frec$Mean

ggplot(over, aes(x=Mean_Ratio, fill=Class)) + 
  geom_histogram(binwidth=0.1)+theme_minimal()+
  xlab("Portion of overshoot sets") + ylab("Number of topologies")+ 
  scale_fill_manual(values=c("DarkMagenta","DarkOrange"))+ theme_minimal() +
  theme(
    panel.grid.major= element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank())

no_overShootTop <- over$Topology[over$Mean_OV==0]

write.csv(no_overShootTop,"no_overShootTop_V2_2.csv")

ggplot(over[over$Mean_BB>=median(over$Mean_BB),], aes(x=Mean_Ratio)) + 
  geom_histogram(binwidth=0.1)+theme_minimal()
ggplot(over[over$Mean_BB<median(over$Mean_BB),], aes(x=Mean_Ratio)) + 
  geom_histogram(binwidth=0.1)+theme_minimal()

ggplot(over, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="fill", stat="identity")




ggplot(over[over$Mean_BB>=mean(over$Mean_BB),], aes(x=Mean_Ratio)) + 
  geom_histogram(binwidth=0.1)+theme_minimal()
ggplot(over[over$Mean_BB<mean(over$Mean_BB),], aes(x=Mean_Ratio)) + 
  geom_histogram(binwidth=0.1)+theme_minimal()

ggplot(over[over$Mean_BB>=30,], aes(x=Mean_Ratio)) +
  geom_histogram(binwidth=0.1)+theme_minimal()
ggplot(over[over$Mean_BB<30,], aes(x=Mean_Ratio)) + 
  geom_histogram(binwidth=0.1)+theme_minimal()


#Analysis of 125 robust topologies

robustTop <- topos[frecb$Topology,]



topNum<-6431
names_vector<- nodes_combi_split_filter2[392,]
test1 <- topos[topNum,]
A <- as.character(names_vector$Node.1[1])
B <- as.character(names_vector$Node.2[1])
C <- as.character(names_vector$Node.3[1])
  namesLinks = c("1:nw -> 1:w","2:n -> 1:e","3:n -> 1:s",
                 "1:ese -> 2:wnw","2:s -> 2:ne","3:e ->2:sw",
                 "1:sw -> 3:nw","2:w -> 3:ne","3:w -> 3:sw")
  positive <- paste(namesLinks[which(test1>0)],collapse="; ")
  negative <- paste(namesLinks[which(test1<0)],collapse="; ")
  null<- paste(namesLinks[which(test1==0)],collapse="; ")
  #png(paste0( "Minimal",toPlot_ID[k],".png"), units="px", width=1500, height=1500, res=600)
  grViz(paste0("
               digraph circles {
               
               # a 'graph' statement
               graph [overlap = true, fontsize = 2, layout = neato]
               
               node [shape = oval,fontsize = 7
               fontname = Helvetica
               fixedsize = true,
               height = 0.2, width = 0.4,peripheries = 1] //sets as circles
               1[label=", A," ","pos=\"0,0.8!\"]
               2[label=", B, " ","pos=\"0.75,0.4!\"]
               3[label=", C, " ","pos=\"0,0!\"]
               
               edge [color = DeepSkyblue2,arrowsize=0.3, arrowhead=normal,  weight =2]",positive, 
               " edge [color = Red4,arrowsize=0.4, arrowhead=tee,weight =2]",negative, 
               "}"))

  names_vector<-combi_names
  topNum <- row.match(temTop,topos)
  topNum
  topNum<-6431
  names_vector<- nodes_combi_split_filter2[h,]
  test1 <- topos[topNum,]
  A <- names_vector[1]
  B <- names_vector[2]
  C <- names_vector[3]
  namesLinks = c("1:nw -> 1:w","2:n -> 1:e","3:n -> 1:s",
                 "1:ese -> 2:wnw","2:s -> 2:ne","3:e ->2:sw",
                 "1:sw -> 3:nw","2:w -> 3:ne","3:w -> 3:sw")
  positive <- paste(namesLinks[which(test1>0)],collapse="; ")
  negative <- paste(namesLinks[which(test1<0)],collapse="; ")
  null<- paste(namesLinks[which(test1==0)],collapse="; ")
  #png(paste0( "Minimal",toPlot_ID[k],".png"), units="px", width=1500, height=1500, res=600)
  grViz(paste0("
               digraph circles {
               
               # a 'graph' statement
               graph [overlap = true, fontsize = 2, layout = neato]
               
               node [shape = oval,fontsize = 7
               fontname = Helvetica
               fixedsize = true,
               height = 0.2, width = 0.4,peripheries = 1] //sets as circles
               1[label=", A," ","pos=\"0,0.8!\"]
               2[label=", B, " ","pos=\"0.75,0.4!\"]
               3[label=", C, " ","pos=\"0,0!\"]
               
               edge [color = DeepSkyblue2,arrowsize=0.3, arrowhead=normal,  weight =2]",positive, 
               " edge [color = Red4,arrowsize=0.4, arrowhead=tee,weight =2]",negative, 
               "}"))
  
  
  row.match(temTop,topos)
  
  
  topNum<-6431
  
  a<- as_data_frame(sub_V1_Net, what = c("edges", "vertices", "both"))
  names_vector<- unique(c(a[,1],a[,2]))
  test1 <- topos[topNum,]
  A <- names_vector[1]
  B <- names_vector[2]
  C <- names_vector[3]
  main_nodes<-paste(combi_names,collapse=" ")
  no_main <- paste(names_vector[!is.element(names_vector,combi_names)],collapse=" ")
  namesLinks = paste0(a[,1],"->",a[,2])
  positive <- paste(namesLinks[which(net_Test$edges$group>0)],collapse="; ")
  negative <- paste(namesLinks[which(net_Test$edges$group<0)],collapse="; ")
  null<- paste(namesLinks[which(test1==0)],collapse="; ")
  #png(paste0( "Minimal",toPlot_ID[k],".png"), units="px", width=1500, height=1500, res=600)
  grViz(paste0("
               digraph circles {
               
               # a 'graph' statement
               graph [overlap = true, fontsize = 8, layout = dot]
               
               node [shape = oval,fontsize = 18
               fontname = Helvetica
               fixedsize = true,
               height = 0.6, width = 1.1,peripheries = 1.5] //sets as circles
               node [fillcolor = lightgrey, style = filled] ",main_nodes,
               " node [fillcolor = white]", no_main,
               " edge [color = DeepSkyblue2,arrowsize=1, arrowhead=normal,penwidth=1.5] ",positive, 
               " edge [color = Red4,arrowsize=1.5, arrowhead=tee,penwidth=1.5] ",negative, 
               "}"))
  
 mermaid(paste0("
               digraph circles {
               
               # a 'graph' statement
               graph [overlap = false, fontsize = 2, layout = neato]
               
               node [shape = oval,fontsize = 7
               fontname = Helvetica
               fixedsize = true,
               height = 0.2, width = 0.4,peripheries = 1] //sets as circles
               
               edge [color = DeepSkyblue2,arrowsize=0.3, arrowhead=normal,  weight =2] ",positive, 
               " edge [color = Red4,arrowsize=0.4, arrowhead=tee,weight =2] ",negative, 
               "}"))
 