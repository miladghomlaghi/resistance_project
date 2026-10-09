library(extrafont)
font_import()
loadfonts(device = "win")
library(readxl)
library(reshape2)
library(ggplot2)
library(dplyr)
load("./Strict/Data/frequency_matlab_210220.RData")
data<-frec


png("Strict/Plots/funtionalVSRebound.png",units="cm", width=9, height=7, res=300)
ggplot(data, aes(x=Mean_GS, y=Mean_BB)) +
  geom_point(color=ifelse(data$Mean>0.03,"#8b008b", "#8b008b"), 
             shape= ifelse(data$Mean>0.01,"circle", "circle"),
             size= ifelse(data$Mean>0.03,1, 1)) +
  #stat_smooth(method = "lm", col = "red")+
  geom_smooth(method = lm, col='black', size=0.6)+
  theme_classic() +
  theme(text=element_text(size=11,colour = "black"),
        panel.border = element_blank(),
        panel.grid = element_blank(), 
        panel.background =  element_rect(fill = "white"),
        #axis.text= element_text(size=10),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.line = element_line(size=0.75, arrow = arrow(angle=18, length =unit(0.15, "inches")))
        # panel.grid.minor.y = element_line(colour="black", size=1),
  ) +xlab("Number of functional sets") +ylab("Number of rebound sets") 
dev.off()

cor(data$Mean_GS,data$Mean_BB)

