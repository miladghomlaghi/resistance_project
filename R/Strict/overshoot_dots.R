
library(readxl)
library(ggplot2)

#data_all<-read_excel("Strict/Data/frequency_replicates_strict.xls",col_names = T)
load("./Strict/Data/frequency_matlab_210220.RData")
names<-c("e11","e21","e31","e12","e22","e32","e13","e23","e33")
data <- over


data$shape<-0
data$shape[data$Mean>0.01 & data$Mean_Ratio>0.20]<-1
data$shape[data$Mean>0.03 | data$Mean_Ratio>0.60 ]<-2
data$shape<-as.factor(data$shape)

data$color<-0
data$color[data$Mean>0.01 & data$Mean_Ratio>0.20]<-1
data$color[data$Mean>0.03]<-2
data$color[data$Mean_Ratio>0.60 ]<-3
data$color<-as.factor(data$color)


png("Strict/Plots/dots_overshoot.png",units="cm", width=18, height=12, res=300)
ggplot(data, aes(x=Mean*100, y=Mean_Ratio*100, shape=shape, color=color)) +
  geom_point(#color=ifelse(data$Mean_Ratio>0.6,"blue", "#8b008b"), 
             #shape=shape,
             #shape= ifelse(data$Mean>0.03,"square", "circle"),
             size= ifelse(is.element(data$shape,c("1","2")),2, 1.5) ) +
  theme_bw()+
  theme(text=element_text(size=12,colour = "black", family="Heveltica"),
        panel.border = element_blank(),
        panel.grid = element_blank(), 
        panel.background =  element_rect(fill = "white"),
        axis.text= element_text(size=12),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.line = element_line(size=1, arrow = arrow(angle=18, length =unit(0.15, "inches")))
       # panel.grid.minor.y = element_line(colour="black", size=1),
  ) +xlab("Rebound score(%)") +ylab("Overshoot score(%)") + 
  scale_color_brewer(palette = "Paired")
dev.off()

overshoot<- data[data$Mean_Ratio>0.6,]

data<-topos[overshoot$Topology,]
names(data)<-names

l<-cluster_topologies(data,2, dist_method="canberra",  hclust_method="complete" )
l$cluster

clust <-row.match(get_clustTop(data,l$groups),topos)
top_Core_Plus <- cbind(clust,100*frecuency_Rep[clust,'Mean'])
top_Core_Plus

png("./Strict/Plots/cluster_overshoot.png", units="cm", width=16, height=16, res=300)
l$cluster
dev.off()


###########
data_IFFL<-topos[frec$Topology,]
data_IFFL_1<- prodlim::row.match(data_IFFL[,c(4,7,8)], t(c(1,1,-1))) 
data_IFFL_2<- prodlim::row.match(data_IFFL[,c(4,7,8)], t(c(-1,1,1))) +1
data_IFFL_1[is.na(data_IFFL_1)] <- 0
data_IFFL_2[is.na(data_IFFL_2)] <- 0
data<-cbind(data, IFFL=as.factor((data_IFFL_1+data_IFFL_2)))

png("Strict/Plots/dots_overshoot_IFFL.png",units="cm", width=18, height=12, res=300)
ggplot(data, aes(x=Mean*100, y=Mean_Ratio*100, shape=shape, color=IFFL)) +
  geom_point(#color=ifelse(data$IFFL==0,"#a6cee3", "#8b008b"), 
    #shape=shape,
    #shape= ifelse(data$Mean>0.03,"square", "circle"),
    size= ifelse(is.element(data$shape,c("1","2")),2, 1.5) ) +
  theme_bw()+
  theme(text=element_text(size=12,colour = "black", family="Heveltica"),
        panel.border = element_blank(),
        panel.grid = element_blank(), 
        panel.background =  element_rect(fill = "white"),
        axis.text= element_text(size=12),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.line = element_line(size=1, arrow = arrow(angle=18, length =unit(0.15, "inches")))
        # panel.grid.minor.y = element_line(colour="black", size=1),
  ) +xlab("Rebound score(%)") +ylab("Overshoot score(%)") + scale_color_viridis(discrete = TRUE, option = "G")
  scale_color_brewer(palette = "Dark2")
dev.off()
                                  