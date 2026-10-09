plot_bar_replicates <-function (frec_Core){

#frec_Core <-frecuency_Rep[core_topos,]

#meltedfrec_Core <- melt(frec_Core[,1:5],id.vars = 1) #topo more than 100, replicates and mean
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
p<-ggplot(my_sum_Core) +
  geom_bar( aes(x=reorder(as.character(Topology),-mean), y=mean),
            stat="identity", fill="DarkOrange3", alpha=0.5, width = 0.4) +
  theme_light() +
  geom_errorbar( aes(x=as.character(Topology), ymin=mean-se, ymax=mean+se),
                 width=0.4, colour="gray26", alpha=0.9, size=0.5) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle=90,size=12,hjust=0.95,vjust=0.25),
    axis.text=element_text(size=18),
    axis.title=element_text(size=18,face="bold"),
    title =element_text(size=16, face='bold')
  ) + 
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 150)) +
  geom_point(data=meltedfrec_Coreb, aes(x=as.character(Topology), y=value), 
             color="black", size=1.5)+
  ylab("Number of sets") + 
  xlab("Topology")

p

#+
# labs(title = "S-sets for core topologies", subtitle = "Mean value of replicates (n=3)")
#+ylim(0, 150)
}
