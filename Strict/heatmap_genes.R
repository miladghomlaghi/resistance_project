library(dplyr)
library(rjson)
Compendium_Cancer_Genes <- read.delim("Compendium_Cancer_Genes.tsv", header=T)

genes_drugs <- unique(Compendium_Cancer_Genes[,c("SYMBOL","CANCER_TYPE", "QVALUE_COMBINATION")])
data<-genes_drugs

json_file <- "drug-ndc-0001-of-0001.json"
json_data <- fromJSON(paste(readLines(json_file), collapse=""))


json_file2 <- "drug-enforcement-0001-of-0001.json"
json_data2 <- fromJSON(paste(readLines(json_file2), collapse=""))


json_file3 <- "drug-event-0040-of-0040.json"
json_data3 <- fromJSON(paste(readLines(json_file3), collapse=""))

min(data$QVALUE_COMBINATION)
max(data$QVALUE_COMBINATION)

initials <- "XYZ"
patern <- paste("^[", initials,"]", sep='')

dataA<- data[grep(patern, data$SYMBOL),]
dataA$SYMBOL<- as.factor(dataA$SYMBOL)
dataA$CANCER_TYPE<- as.factor(dataA$CANCER_TYPE)
dataA$SYMBOL<- factor(dataA$SYMBOL, levels=rev(levels(dataA$SYMBOL)))

# all possible combinations
#all <- dataA %>% expand(SYMBOL, CANCER_TYPE)

# join with all, n will be NA for obs. in all that are not present in v
#dataA = dataA%>% group_by_at(vars(SYMBOL, CANCER_TYPE)) %>% 
 # summarize(QVALUE_COMBINATION = mean()) %>% right_join(all)


file_name = paste("Strict/Plots/genes_cancer",initials,'.png', sep='')
png(file_name,units="cm", width=14, height=15, res=300)

ggplot(dataA, aes(CANCER_TYPE, SYMBOL, fill= QVALUE_COMBINATION)) + 
  geom_tile(color="black",
            size=0.25) +
 # scale_fill_distiller(palette = "RdPu",trans = "reverse") +
#  theme_ipsum()+
  theme_bw(base_size=8)+
  theme(text=element_text(size=9,colour = "black"),
        axis.ticks=element_line(size=0.4),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.placement = "outside",
        plot.background=element_blank(),
        legend.position = "top",
        #remove plot border
      #  axis.title.y = element_blank(),
      #  axis.title.x = element_blank(),
        panel.grid.major = element_blank())+ xlab("Cancer type") + ylab("Driver gene")+
        #strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
  #ggtitle(label = "Microbe Class Abundance") + xlabr\

  scale_fill_gradientn(name = "Mean q-value",
                      limits = c(min(data$QVALUE_COMBINATION),0.5),
                     breaks = c(min(data$QVALUE_COMBINATION),0.05,0.5),
                      labels = c("","0.05" ,"0.5"),
                      colours = c("#006DAE","#92b9d4","#a7cee8","#d4ebfa","#e3f2fc","#f5fbff",
                     "#FFFFFF"))
    
dev.off()


ggplot(m3,aes(x=year,y=state,fill=count))+
  #add border white colour of line thickness 0.25
  geom_tile(colour="white",size=0.25)+
  #remove x and y axis labels
  labs(x="",y="")+
  #remove extra space
  scale_y_discrete(expand=c(0,0))+
  #define new breaks on x-axis
  scale_x_discrete(expand=c(0,0),
                   breaks=c("1930","1940","1950","1960","1970","1980","1990","2000"))+
  #set a base size for all fonts
  theme_grey(base_size=8)+
  #theme options
  theme(
    #bold font for legend text
    legend.text=element_text(face="bold"),
    #set thickness of axis ticks
    axis.ticks=element_line(size=0.4),
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank())


###########
#thesis
ggplot(dataA, aes(CANCER_TYPE, SYMBOL, fill= value)) + 
  geom_tile(colour="black") +
  # scale_fill_distiller(palette = "RdPu",trans = "reverse") +
  theme_grey(base_size=8)+
  theme(text=element_text(size=7.5,colour = "black"),
        axis.ticks=element_line(size=0.4),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.placement = "outside",
        plot.background=element_blank(),
        #remove plot border
        panel.border=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  #strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
  #ggtitle(label = "Microbe Class Abundance") +
  
  scale_fill_gradient(name = "Mean q-value",
                      low = "#006DAE",
                      high = "#F6F6F6", trans = 'log')

dev.off()
ggplot(dataA, aes(CANCER_TYPE, SYMBOL, fill= value)) + 
  geom_tile(colour="black") +
  # scale_fill_distiller(palette = "RdPu",trans = "reverse") +
  theme_grey(base_size=8)+
  theme(text=element_text(size=7.5,colour = "black"),
        axis.ticks=element_line(size=0.4),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.placement = "outside",
        plot.background=element_blank(),
        #remove plot border
        panel.border=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  #strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
  #ggtitle(label = "Microbe Class Abundance") +
  
  scale_fill_gradient(name = "Mean q-value",
                      low = "#006DAE",
                      high = "#F6F6F6", trans = 'log')

dev.off()


png("Strict/Plots/genes_cancerA.png",units="cm", width=14, height=18, res=300)

ggplot(dataA, aes(Var2, Var1, fill= Freq)) + 
  geom_tile(colour="white",size=0.25) +
  # scale_fill_distiller(palette = "RdPu",trans = "reverse") +
  #  theme_ipsum()+
  theme_grey(base_size=8)+
  theme(text=element_text(size=9,colour = "black"),
        axis.ticks=element_line(size=0.4),
        axis.text.x = element_text(size=7,angle = 90, vjust = 0.5, hjust=1),
        strip.placement = "outside",legend.position = "none",
        plot.background=element_blank(),
        #remove plot border
        panel.border=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  #strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
  #ggtitle(label = "Microbe Class Abundance") +
  
  scale_fill_gradient(name = "Sqrt(Abundance)",
                      low = "#FFFFFF",
                      high = "#006DAE")

dev.off()






