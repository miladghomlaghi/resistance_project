library(extrafont)
font_import()
loadfonts(device = "win")
library(readxl)
library(reshape2)
library(ggplot2)
library(dplyr)
load("./Strict/Data/frequency_matlab_210220.RData")
data<-frecuency_Rep[frecuency_Rep$Mean>0.01, ]


# frec <-frequency_Rep[!frequency_Rep$Max==0,]
# frecb <-frequency_Rep[frequency_Rep$Mean>100,]


meltedFreb <- melt(data[,c("Topology", "Rep_1", "Rep_2","Rep_3")],id.vars = 1) #topo more than 100, replicates and mean
meltedFreb$value<-meltedFreb$value*100

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
png("./Strict/Plots/BarDeviation_thesis.png", units="cm", width=18, height=9, res=300)
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
    axis.text.y = element_text(size=10),
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

