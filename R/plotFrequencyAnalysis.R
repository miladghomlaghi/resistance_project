R

install.packages('numbers')
install.packages('gstat')


# Library
library(tidyverse)
#library(ggridges)
library(ggplot2)
library(reshape2)
library(scales)
library(gridExtra)
library(RColorBrewer)
library(rstudioapi)
library(reshape2)
library(ggplot2)
#library(ggridges)
#library(geosphere)
library(gstat)
library(windfarmGA)
library(viridis)
library(numbers)
library(tidyverse)
library(foreach)
#library(BBmisc)
library(RColorBrewer)

# frequency_Rep <- read.csv("C:/Users/kisl1/Downloads/PhD/RCodes/frecuency_Rep.csv", header=FALSE)
# foreach (j = 1:nrow(frequency_Rep)) %do%  (frequency_Rep[j,5]<-mean(as.numeric(frequency_Rep[j,2:4])))
# 
# foreach (j = 1:nrow(frequency_Rep)) %do%  (frequency_Rep[j,6]<-max(as.numeric(frequency_Rep[j,2:4])))
# 
# frequency_Rep <- cbind(frequency59040_R1Seq,frequency59040_R1Seq[,2],frequency59040_R1Seq[,2])
# frequency_Rep<-cbind(frequency_Rep, rowMeans(frequency_Rep[,2:4]))
# 
# 
# names(frequency_Rep) <- c("Topology", "Rep_1", "Rep_2","Rep_3", "Mean")
# frec <-frequency_Rep[!frecuency_Rep$Max==0,]

frecb <-frequency_Rep[frequency_Rep$Mean>0,]
frec <-frecb


meltedFre <- melt(frec[,1:5],id.vars = 1) #all topologies
meltedFreb <- melt(frecb[,1:5],id.vars = 1) #topo more than 100, replicates and mean
meltedFre2b <- melt(frecb[,1:4],id.vars = 1) #topo more thatn 100, only mean


data=meltedFre[meltedFre$variable=="Mean",]

##### Sorted

data <-data[order(-data$value),]

png("sorted1.png", units="px", width=1805, height=1200, res=300)
ggplot(data, aes(reorder(Topology, -value), value)) +
  geom_segment(aes(x=reorder(Topology, -value), xend=reorder(Topology, -value),y=0, yend=value),color=ifelse(data$value > 100,"orangered3","darkgrey")) + 
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=16),
    axis.title=element_text(size=16,face="bold"),
    axis.title.x = element_text(margin = margin(t = -20)),
    title =element_text(size=16, face='bold')
  ) +
  ylab("Number of sets") + xlab("Sorted topologies") +
  labs(title = "S-sets for each topology", subtitle = "Mean value of three replicates")+
  annotate("text", x = 200, y = 110, label = "1% of 10,000 sets", color="orangered3", size=6 , angle=0, fontface="bold",  hjust=0)

# insert ggplot code
dev.off()


##### No Sorted
png("NoSorted1.png", units="px", width=1805, height=1200, res=300)

ggplot(data, aes(x=Topology, y=value)) +
  geom_segment(aes(x= Topology, xend=Topology,y=0, yend=value),color="darkgrey") +
  geom_point(size=ifelse(data$value > 5, 1.3, 1), color=ifelse(data$value > 5,"orangered3","dodgerblue")) + theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text=element_text(size=18),
    axis.title=element_text(size=18,face="bold"),
    title =element_text(size=20, face='bold')
  ) +
  xlab("Topology #id") +
  ylab("Number of sets") + 
  labs(title = "S-sets for each topology", subtitle = "Mean value of replicates (n=3)") +
  xlim(0, 59100)

dev.off()

p <- data %>%
  filter( price<300 ) 
  


tmp_size<-nrow(data)
for (i in 1:10){
data[(tmp_size+1+(i-1)*50):(tmp_size+50*i),3]<-i
data[(tmp_size+1+(i-1)*50):(tmp_size+50*i),1]<-0
}
data[,2]<-"Mean"
data$value<-as.factor(data$value)


data$value

  ggplot(data, aes(x=value)) +
  geom_histogram( stat="count", fill="#69b3a2", color="#e9ecef", alpha=0.9) +
    theme_classic()+
    stat_count(geom='text',color='black', aes(label=(..count..)-50), size=5,
               position=position_dodge(width=0.9), vjust=-0.25)  +  ylab("Number of topologies")+
  theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(size=18),
      axis.text=element_text(size=18),
      axis.title=element_text(size=18,face="bold"),
      title =element_text(size=18, face='bold'),
      text = element_text(size = 18)) +scale_y_continuous(expand = c(0,0),limits=c(0,10000)) + 
    xlab("Number of sets") 
  
  
  
  theme_ipsum() +

# Calculates mean, sd, se and IC
my_sum <- meltedFre2b%>%
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
  geom_bar( aes(x=reorder(as.character(Topology),-mean), y=mean), stat="identity", fill="firebrick", alpha=0.5, width = 0.6) +
  theme_light() +
  geom_errorbar( aes(x=as.character(Topology), ymin=mean-sd, ymax=mean+sd), width=0.4, colour="gray26", alpha=0.9, size=0.5) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle=90,size=10,hjust=0.95,vjust=0.25),
    axis.text=element_text(size=18),
    axis.title=element_text(size=18,face="bold"),
    title =element_text(size=18, face='bold')
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 175)) +
  geom_point(data=meltedFre2b, aes(x=as.character(Topology), y=value), 
             color="dodgerblue4", size=1)+
  ylab("Number of sets") + 
  xlab("Topology")+
  labs(title = "S-sets for robust topologies", subtitle = "Three replicates") 



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
  labs(title = "S-sets for core topologies", subtitle = "Mean value of replicates (n=3)") 
