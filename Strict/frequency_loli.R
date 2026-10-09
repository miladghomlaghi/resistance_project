
library(readxl)

data_all<-load("./Strict/Data/frequency_matlab_210220.RData")
data<-frecuency_Rep[frecuency_Rep$Mean>0, ]
#data<- melt(data[,c("ID", "Rep_1", "Rep_2","Rep_3", "Mean")],id.vars = 1)

##### No Sorted
png("./Strict/Plots/freq_loli.png", units="cm", width=24, height=16, res=300)
ggplot(data, aes(x=Topology, y=Mean*100)) +
  geom_segment(aes(x= Topology, xend=Topology,y=0, yend=Mean*100),color=ifelse(data$Mean > 0.5,"#8b008b","#8b008b"))+
  geom_point(size=ifelse(data$Mean > 0.25, 2, 2), color=ifelse(data$Mean > 0.5,"#8b008b","#8b008b")) +
  theme(text=element_text(size=20,colour = "black", family="Heveltica"),
        #panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        panel.grid = element_blank(), 
        panel.background =  element_rect(fill = "white"),
        axis.text.y = element_text(size=20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_line(size=1, arrow = arrow(angle=18, length =unit(0.15, "inches"))),
       panel.grid.minor.y = element_line(colour="black", size=1),
  )+
  xlab("Topologies per ID") +
  ylab("Rebound score (%)") + scale_y_continuous(minor_breaks = c(1,3),limits=c(0,6)) +xlim(0, 16038) 

  #labs(title = "S-sets for each topology", subtitle = "Mean value of replicates (n=3)")
dev.off()
