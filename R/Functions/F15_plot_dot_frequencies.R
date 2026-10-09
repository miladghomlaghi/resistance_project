plot_dot_frequencies<-function (frequency_Rep, given_topology){

 
  frec <-frequency_Rep[frequency_Rep$Mean>0,]
  meltedFre <- melt(frec[,1:5],id.vars = 1) #all topologies
  data=meltedFre[meltedFre$variable=="Mean",]  
  

if(identical(which(data$Topology == given_topology),integer(0))){
  data<-rbind(data,data[nrow(data),])
  data[nrow(data),1]<-given_topology
  data[nrow(data),3]<-0
}
p<-ggplot(data, aes(x=Topology, y=value)) +
  geom_segment(aes(x= Topology, xend=Topology,y=0, yend=value),color="#7a7d8c", size=1)+
  geom_point(size=ifelse(data$Topology == given_topology, 2.5, 1), 
             color=ifelse(data$Topology == given_topology,'#e1af28',"#a503ab")) + theme_classic() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text=element_text(size=10),
    axis.title=element_text(size=10,face="bold"),
    title =element_text(size=10, face='bold')
  ) +
  xlab("Network ID") +
  ylab("Score") + 
  labs(title = "") +
  xlim(0, nrow(frequency_Rep)) 
return(p)
}
