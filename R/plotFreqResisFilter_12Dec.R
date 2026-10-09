

# Library
library(tidyverse)
library(ggridges)
library(ggplot2)
library(reshape2)
library(scales)
library(gridExtra)
library(RColorBrewer)
library(rstudioapi)
library(reshape2)
library(ggridges)
library(geosphere)
library(gstat)
library(windfarmGA)
library(viridis)
library(numbers)
library(tidyverse)
library(foreach)
library(BBmisc)
library(dplyr)

# frequency_Rep <- read.csv("C:/Users/kisl1/Downloads/PhD/RCodes/frecuency_Rep.csv", header=FALSE)
# foreach (j = 1:nrow(frequency_Rep)) %do%  (frequency_Rep[j,5]<-mean(as.numeric(frequency_Rep[j,2:4])))
# 
# # foreach (j = 1:nrow(frequency_Rep)) %do%  (frequency_Rep[j,6]<-max(as.numeric(frequency_Rep[j,2:4])))
# 
# 
# names(frequency_Rep) <- c("Topology", "Rep_1", "Rep_2","Rep_3", "Mean","Max")
# frec <-frequency_Rep[!frequency_Rep$Max==0,]
# frecb <-frequency_Rep[frequency_Rep$Mean>100,]


meltedFre <- melt(frec[,1:5],id.vars = 1) #all topologies
meltedFreb <- melt(frecb[,c("Topology", "Rep_1", "Rep_2","Rep_3", "Mean")],id.vars = 1) #topo more than 100, replicates and mean
meltedFre2b <- melt(frecb[,1:4],id.vars = 1) #topo more thatn 100, only mean


data=meltedFre[meltedFre$variable=="Mean",]






##### No Sorted
png("NoSorted1.png", units="cm", width=24, height=16, res=300)
ggplot(data, aes(x=Topology, y=value)) +
  geom_segment(aes(x= Topology, xend=Topology,y=0, yend=value),color=ifelse(data$value > 3,"#8b008b","#8b008b"))+
  geom_point(size=ifelse(data$value > 3, 2, 1.5), color=ifelse(data$value > 3,"#8b008b","#8b008b")) +
  theme(text=element_text(size=12,colour = "black", family="Heveltica"),
        #panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        panel.grid = element_blank(), 
        panel.background =  element_rect(fill = "white"),
        axis.text.y = element_text(size=12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_line(size=1, arrow = arrow(angle=18, length =unit(0.15, "inches")))
  )+
  xlab("Topology per ID") +
  ylab("Performance (%)") + 
  #labs(title = "S-sets for each topology", subtitle = "Mean value of replicates (n=3)") +
  xlim(0, 16038)

dev.off()


##### No Sorted, No text
#png("NoSorted1.png", units="px", width=1805, height=1200, res=300)

ggplot(data, aes(x=Topology, y=value)) +
  geom_segment(aes(x= Topology, xend=Topology,y=0, yend=value),color="darkgrey", size=1)+
  geom_point(size=ifelse(data$value > 100, 1.5, 1), color=ifelse(data$value > 100,"orangered3","dodgerblue")) + theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text=element_text(size=18),
    axis.title=element_blank()
    #title =element_text(size=20, face='bold')
  ) + xlim(0, 16038)
#dev.off()

frecGS<-frequency59040_R1Seq
names(frecGS)<- c('Topology','Mean_GS')
meltedFrec <- melt(frecGS)

##### Sorted

#data <-data[order(-data$value),]


##### Scatter plot functional sets for each topology
png("ScatterFuntional.png", units="cm", width=24, height=16, res=300)
ggplot(frecGS, aes(x=Topology, y=Mean_GS)) +
  geom_point(size=0.6, 
             color=ifelse(frecGS$Mean_BB == 0,"SteelBlue","DarkMagenta")) + theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text=element_text(size=12),
    axis.title=element_text(size=12,face="bold"),
    title =element_text(size=12, face='bold')
  ) +
  xlab("Topology ID") +
  ylab("") + scale_y_continuous(limits=c(0, 5000), position = "right")+
  labs(title = "Functional parameters per topology") +
  xlim(0, 16038)
dev.off()


##### Scatter plot Relation between functional sets and bb-sets
png("FunctionalVSbb.png", units="px", width=1800, height=1200, res=300)
ggplot(frecGS, aes(x=Mean_GS,y=Mean_BB)) +
  geom_point(size=0.5, 
             color=ifelse(frecGS$Mean_BB == 0,"SteelBlue","DarkMagenta")) + theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text=element_text(size=12),
    axis.title=element_text(size=12),
    title =element_text(size=12)
  ) +
  xlab("Functional sets") +
  ylab("Bounce-back sets") + scale_y_continuous(limits=c(0, 120), position = "right")+
  labs(title = "") +
  xlim(0, 4500)
dev.off()


##### Scatter plot Relation between functional sets and bb-sets
png("bbVSFunctional.png", units="px", width=1800, height=1200, res=300)
ggplot(frecGS, aes(x=Mean_BB,y=Mean_GS)) +
  geom_point(size=0.5, 
             color=ifelse(frecGS$Mean_BB == 0,"SteelBlue","DarkMagenta")) + theme_classic() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text=element_text(size=11),
    axis.title=element_text(size=12),
    title =element_text(size=12)
  ) +
  xlab("Bounce-back sets") +
  ylab("Functional sets") + scale_y_continuous(limits=c(0, 4500), position = "left")+
  labs(title = "") +
  xlim(0, 125)
dev.off()



# Desinty of funtional sets
frecGS$Type <- ifelse(frecGS$Mean_BB==0,"A","B")
png("densitySets.png", units="px", width=1200, height=1200, res=300)
ggplot(frecGS, aes(x=Mean_GS, color=Type, fill=Type)) +
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5,show.legend = FALSE)+
  geom_density(alpha=0.6,show.legend = FALSE)+
  scale_color_manual(values=c("DarkMagenta", "SteelBlue", "DarkMagenta"))+
  scale_fill_manual(values=c("DarkMagenta", "SteelBlue", "DarkMagenta"))+
  labs(title="",x="Funtional parameter sets", y = "Density")+ 
  theme_classic()
dev.off()



#png("sorted1.png", units="px", width=1805, height=1200, res=300)
ggplot(meltedFre, aes(reorder(Topology, -value), value)) +
  geom_segment(aes(x=reorder(Topology, -value), xend=reorder(Topology, -value),y=0, yend=value),
               color="GainsBoro",size=0.2)+
  geom_point(size=ifelse(meltedFre$value < 1, 1.3, 2), 
             color=ifelse(meltedFre$value < 1,"DarkMagenta","DarkOrange")) + 
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=20),
    axis.title=element_blank(),
    axis.title.x = element_blank(),
    title =element_blank()
  ) +
  ylab("") + xlab("") +ylim("1","2","3","4","5")+
  labs(title = "", subtitle = "")+ scale_x_discrete(expand=expand_scale(0.02))

+
  annotate("text", x = 400, y = 80, label = "At least 1% of 5000 sets", 
           color="DarkOrange", size=6 , angle=0, hjust=0) 

top_Core <- cbind(coreTopIDs,frecuency_Rep[coreTopIDs,'Mean'])
good_clus <-  c(1,2,3,5,6,8,11,12,14,18,20,24,27,28,31,32,33,36,44,45,46,47,48)
top_Core <- top_Core[good_clus,]
top_Core$Mean <-top_Core$Mean*100
clus_NoSets <- as.numeric(rownames(top_Core[top_Core$Mean==0,]))
clus_NoSets

good_clus <- good_clus[-match(clus_NoSets,good_clus)]


##### No Sorted dodgerblue

data_highlight <- results_UU
data_highlight <- cbind(coreTopIDs_Plus,frecuency_Rep[coreTopIDs_Plus,"Mean"]*100)
data_highlight <- cbind(c(5464,6434,8213,8039),frecuency_Rep[c(5464,6434,8213,8039),"Mean"]*100) #two iterations for clustering of plus V2_2
names(data_highlight)[1]<-"Topology"
"DarkOrange"

png("FrequencyPercentaje.png", units="px", width=1800, height=1200, res=300)
ggplot(meltedFre, aes(x=Topology, y=value)) +
  geom_segment(aes(x= Topology, xend=Topology,y=0, yend=value),color="darkgrey",size=0.3)+
  geom_point(size=ifelse(meltedFre$value > 0.99, 0.8, 0.5), 
             color=ifelse(meltedFre$value < 3,"DarkMagenta","DarkOrange")) + theme_classic() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text=element_text(size=10),
    axis.title=element_text(size=10,face="bold"),
    title =element_text(size=10, face='bold')
  ) +
  xlab("Topology ID") +
  ylab("Parameter sets (%)") + 
  labs(title = "") +
  xlim(0, 16038) 
dev.off()

+ geom_point(data = data_highlight, aes(x=Topology, y=Mean,color=Degree), size=1) 
                              

#+ geom_point(data = data_highlight, aes(x=Topology, y=Mean), size=1, 
 #            color='Chartreuse3') 

+
  geom_point(aes(x=Topology, y=frecGS$Mean_BB),size=ifelse(frecGS$Mean > 0, 1, 1.2), 
             color=ifelse(frecGS$Mean_GS > 0,"DarkMagenta","DodgerBlue")) 


geom_segment(aes(x= Topology, xend=Topology,y=0, yend=Mean_BB),color="DarkMagenta")+
  geom_segment(aes(x= Topology, xend=Topology,y=frecGS$Mean_BB, yend=Mean_GS),color="DarkOrange")

# Calculates mean, sd, se and IC
#  newRobTop<- meltedFreb[meltedFreb$value>1 & meltedFreb$variable=="Mean",]
# newRobTop<-newRobTop$Topology
meltedFreb$value <- meltedFreb$value *100
  #meltedFreb[is.element(meltedFreb$Topology,newRobTop),]
my_sum <- meltedFreb%>%
  group_by(Topology) %>%
  summarise( 
    n=n(),
    mean=mean(value),
    sd=sd(value)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

# Standard deviation
png("BarDeviation.png", units="cm", width=12, height=8, res=300)
ggplot(my_sum) +
  geom_bar( aes(x=reorder(as.character(Topology),-mean), y=mean), stat="identity", 
            fill=ifelse(my_sum$mean > 3,"#70ad47","#efc000"),width = 0.5,alpha=1) +
  theme_classic() +
  geom_errorbar( aes(x=as.character(Topology), ymin=mean-sd, ymax=mean+sd), width=0.4, colour="gray26", alpha=0.9, size=0.6)+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=14),
    axis.title.x=element_blank(),
    axis.title=element_text(size=12,face="bold"),
    title =element_text(size=12, face='bold')
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 6)) +
  geom_point(data=meltedFreb, aes(x=as.character(Topology), y=value), 
             color="black", size=0.8)+
  ylab("") + 
  xlab("")+
  labs(title = "") 
dev.off()

#element_text(angle=90,size=16,hjust=0.95,vjust=0.25)

##### Topos ID more that 100 (mean of replicates)

toposID = t(my_sum[,1])




##Core Topos

core_topos = c(701,	1997,	8231,	8960,	9284,	8636,	11471,	17462,	20378,	26453)
frec_Core <-frecuency_Rep[core_topos,]

meltedfrec_Core <- melt(frec_Core[,1:5],id.vars = 1) #topo more than 100, replicates and mean
meltedfrec_Coreb <- melt(frec_Core[,1:4],id.vars = 1) 

# Calculates mean, sd, se and IC
my_sum_Core <- meltedfrec_Coreb%>%
  group_by(Topology) %>%
  summarise( 
    n=n(),
    mean=mean(value),
    sd=sd(value)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

# Standard deviation
ggplot(my_sum_Core) +
  geom_bar( aes(x=reorder(as.character(Topology),-mean), y=mean), stat="identity", fill="firebrick", alpha=0.5, width = 0.4) +
  theme_light() +
  geom_errorbar( aes(x=as.character(Topology), ymin=mean-sd, ymax=mean+sd), width=0.4, colour="gray26", alpha=0.9, size=0.5) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle=90,size=12,hjust=0.95,vjust=0.25),
    axis.text=element_text(size=18),
    axis.title=element_text(size=18,face="bold"),
    title =element_text(size=16, face='bold')
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 150)) +
  geom_point(data=meltedfrec_Coreb, aes(x=as.character(Topology), y=value), 
             color="dodgerblue4", size=1)+
  ylab("Number of sets") + 
  xlab("Topology")+
  labs(title = "S-sets for core topologies", subtitle = "Mean value of replicates (n=3)") +ylim(0, 150)


#png("sorted1.png", units="px", width=1805, height=1200, res=300)
# ggplot(data, aes(reorder(Topology, -value), value)) +
#   geom_segment(aes(x=reorder(Topology, -value), xend=reorder(Topology, -value),y=0, yend=value),
#                color="GainsBoro")+
#   geom_point(size=ifelse(data$value > 50, 1.25, 1), 
#              color=ifelse(data$value < 50,"DarkMagenta","DarkOrange")) + 
#   
#   theme_light() +
#   theme(
#     panel.grid.major.x = element_blank(),
#     panel.border = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.text.x=element_blank(),
#     axis.text.y=element_text(size=16),
#     axis.title=element_text(size=16,face="bold"),
#     axis.title.x = element_text(margin = margin(t = -20)),
#     title =element_text(size=16, face='bold')
#   ) +
#   ylab("") + xlab("") +
#   labs(title = "", subtitle = "")+
#   annotate("text", x = 100, y = 80, label = "At least 1% of 5000 sets", 
#            color="DarkOrange", size=6 , angle=0, hjust=0) 

# insert ggplot code
#dev.off() 

##### overshoot
png("NoSorted1.png", units="px", width=1805, height=1200, res=300)
ggplot(over, aes(x=Topology, y=Mean_OV)) +
  geom_point(color=ifelse(over$Mean_BB < 25,"CornflowerBlue","BlueViolet"), 
             size=ifelse(over$Mean_BB < 25,0.9,1.2))+ 
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    
    axis.text=element_text(size=16),
    axis.title=element_text(size=16,face="bold"),
    title =element_text(size=16, face='bold')
  ) +
  xlab("Topology ID") + ylab("Number of sets")


 + scale_y_continuous(limits=c(1000, 5500), position = "right")+
  labs(title = "Parameter sets per topology that satisfy i") +
  xlim(0, 16038)
dev.off() 


topos<-topos_all[frec[,1],]
namesLinks = c("A-A","B-A","C-A","A-B","B-B","C-B","A-C","B-C","C-C")
namesLinks = c("1-1","2-1","3-1","4-1","1-2","2-2","3-2","4-2","1-3","2-3","3-3","4-3","1-4","2-4","3-4","4-4" )
namesLinks=as.character(namesLinks)
colnames(topos) = namesLinks

robustTop<-frec[,1]
toposMap <-topos_all[unlist(robustTop),] 
rownames(toposMap) <- unlist(robustTop)
colnames(toposMap) = namesLinks

toposMap$Topology<-rownames(toposMap)
toposMelt <-melt(toposMap)
names(toposMelt)[2:3]<-c("Regulation","Value") 

dend.order <- order.dendrogram(dend)
#toposMap2<-toposMap
topos_data2<-topos_data*perce
topos_data2$Topology<-rownames(topos_data)
if(double==1){
  topos_data2[,10:18]<-(-1)*topos_data2[,10:18]}
levelsT = topos_data2$Topology[dend.order]
topos_data2<-topos_data2[dend.order,]
toposMelt <-melt(topos_data2)
names(toposMelt)[2:3]<-c("Regulation","Value") 

# Create dendrogram plot

#toposMelt$Topology <- factor(toposMelt$Topology, levels = levelsT,ordered = TRUE)
toposMelt$Value <- as.factor(toposMelt$Value)
colorLH<-c('Red4','DeepSkyblue2','White')
#colorLH<-c('FireBrick','DodgerBlue')
colorM<-list(color=brewer.pal(9, "Pastel1"))
colorM<-colorM$color[9]

png("HeatNoClus.png", units="px", width=1200, height=1200, res=300)
ggplot(data = toposMelt, aes(x = Regulation, y = Topology)) +
  geom_tile(aes(fill = Value))+
  theme_linedraw(base_size = 14) + scale_y_discrete(position = "left")+
  theme(legend.position="none",
        axis.title.y=element_text(angle=90,vjust=1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = "Black",
                                        colour = "Black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "Black"))



scale_fill_gradient2(name = 'Strength',low=colorLH[1],high=colorLH[2], 
                     mid=colorLH[3],breaks=c(-1,0,1),labels=c("Inhibition",
                                                              "No connection","Activation"))

,
                     limits=c(-1,1)) 


dev.off()
+
  removeGrid()

+ theme_linedraw()

ggplot(data = toposMelt, aes(x = Regulation, y = Topology)) +
geom_rect(aes(fill = Value), colour = "grey50")


################################
###############################

frequency_Rep <- cbind(frequency59040_R1Seq, frequency59040_R1Seq[,2], frequency59040_R1Seq[,2])
names(frequency_Rep)<- c('Topology','Rep_1','Rep_2','Rep_3')

frequency_Rep$Mean <- Means(frequency_Rep[,2:4])
meltedfrec <- melt(frequency_Rep[,1:4],id.vars = 1) 

# Calculates mean, sd, se and IC
my_sum <- meltedfrec%>%
  group_by(Topology) %>%
  summarise( 
    n=n(),
    mean=mean(value),
    sd=sd(value)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

# Standard deviation
ggplot(my_sum) +
  geom_bar( aes(x=reorder(as.character(Topology),-mean), y=mean), stat="identity", fill="firebrick", alpha=0.5, width = 0.4) +
  theme_light() +
  geom_errorbar( aes(x=as.character(Topology), ymin=mean-sd, ymax=mean+sd), width=0.4, colour="gray26", alpha=0.9, size=0.5) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text=element_text(size=18),
    axis.title=element_text(size=18,face="bold"),
    title =element_text(size=16, face='bold')
  ) 

+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 150)) +
  geom_point(data=meltedfrec, aes(x=as.character(Topology), y=value), 
             color="dodgerblue4", size=1)+
  ylab("Number of sets") + 
  xlab("Topology")+
  labs(title = "S-sets for core topologies", subtitle = "Mean value of replicates (n=3)") +ylim(0, 150)




