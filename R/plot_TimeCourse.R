
### Plot time course

x3values<- x3values_2240_6832R1
x4values<- x4values_2240_6832R1
top_num<-2240
set_num<-6832
nodes_num<-4

names(x3values) <- c('D1+D2','D1->D2','D2->D1','D1','D2')
names(x4values) <- c('D1+D2','D1->D2','D2->D1','D1','D2')

x3values<-x3values[,c(4,5,1,2,3)]
x4values<-x4values[,c(4,5,1,2,3)]

#x3values$Time <-  seq(1,60, 0.118)[-1]
x3values$Time <-  seq(1,60)
x3values<-melt(x3values, id.vars = "Time")

#x4values$Time <- seq(1,60, 0.118)[-1]
x4values$Time <- seq(1,60)
x4values<-melt(x4values, id.vars = "Time")


A<- ggplot(x3values,aes(x=Time, y=value, group=variable, color= variable)) +
  geom_line(aes(linetype=variable), size=0.5)+theme_bw()+
  scale_linetype_manual(values=c("dashed","dashed","solid","solid","solid"))+
  scale_color_manual(values=c(rgb(0,153/255,0/255),rgb(0,0/255,255/255),"black",
                              rgb(255/255,170/255,0/255),rgb(144/255,57/255,229/255)))+
  ylim(c(0,1))+ylab("x3 concentration") + scale_fill_discrete(name = "New Legend Title")+
  geom_line(x3values[x3values$variable == "D1+D2",], 
             mapping = aes(x=Time, y=value), size=0.5)

legend <- A
A <- A + theme(legend.position="none")

B<-ggplot(x4values,aes(x=Time, y=value, group=variable, color= variable)) +
  geom_line(aes(linetype=variable), size=0.5)+theme_bw()+
  theme(
        legend.position="none")+
  scale_linetype_manual(values=c("dashed","dashed","solid","solid","solid"))+
  scale_color_manual(values=c(rgb(0,153/255,0/255),rgb(0,0/255,255/255),"black",
                              rgb(255/255,170/255,0/255),rgb(144/255,57/255,229/255)))+
  #scale_shape_manual(values=c(20, "none"))+
  ylim(c(0,1))+ylab("x4 concentration") + 
  geom_line(x4values[x4values$variable == "D1+D2",], 
             mapping = aes(x=Time, y=value), size=0.5)


png(paste0('./Images/',nodes_num,'Nodes/timeCourseX3_T',
           top_num,'_S',set_num,'.png',collapse = ''),
    units="px",width=2443,height = 1651,res=600)
A
dev.off()


png(paste0('./Images/',nodes_num,'Nodes/timeCourseX4_T',
           top_num,'_S',set_num,'.png',collapse = ''),
    units="px",width=2443,height = 1651,res=600)
B
dev.off()



png(paste0('./Images/',nodes_num,'Nodes/legend',
           top_num,'_S',set_num,'.png',collapse = ''),
    units="px",width=2443,height = 1651,res=600)
legend
dev.off()





prow <- plot_grid( A, B, legend, labels = c("A", "B", ""),hjust = -0.5,
                   nrow = 2, rel_widths = c(2.3, 2.3, 0.8))
prow
+ theme(legend.position="none"),
                  p2 + theme(legend.position="none"),
                  p3 + theme(legend.position="none"),
                  align = 'vh',
                  labels = c("A", "B", "C"),
                  hjust = -1,
                  nrow = 1
)

grid.arrange( grobs=lapply(list(A, B, legend), grobTree))
,ncol=3,widths=c(2.3, 2.3, 0.8))

dev.off() 

test<-lapply(list(A, B, legend), grobTree)


theme_classic()+
windowsFonts()
theme_classic()+

142, 0, 229


  +ash", "dotted","longdash")) 
  
  +
  scale_color_manual(values = c("black","#69b3a2
 scale_linetype_manual(values=c("dashed","solid","twod
    scale_size_manual(values=c(0.7,1.3)) ")) +

x3values%>% 
  mutate( highlight=ifelse(variable == "D1->D2", "Sequential","Other")) %>%
  
  
  ggplot(data=x3values, aes(x=Time, y=value, group=variable)) + geom_point()
(aes(type=variable)) 

+ theme_classic()+
  scale_linetype_manual(values=c("dashed","solid","twodash", "dotted","longdash")) 


+
  scale_color_manual(values = c("black","#69b3a2")) +
  scale_size_manual(values=c(0.7,1.3)) 




+
   
+
  

   scale_color_manual(values = c("#00AFBB", "#E7B800","#69b3a2", "#FC4E07","darkorchid1"))+
   
    
    
  +
    scale_linetype_discrete(legend=F) 
  
  
     
    theme(legend.position="none") +
    ggtitle("Popularity of American names in the previous 30 years") +
    +
    geom_label( x=1990, y=55000, label="Amanda reached 3550\nbabies in 1970", size=4, color="#69b3a2") +
    theme(
      legend.position="none",
      plot.title = element_text(size=14)
  +
  scale_fill_viridis(discrete = TRUE) +
  theme(legend.position="none") +
  ggtitle("Popularity of American names in the previous 30 years") +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.stenpacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8),
    plot.title = element_text(size=14)
  ) +
  facet_wrap(~name)