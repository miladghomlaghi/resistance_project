#Figure 1, introduction

library(extrafont)
font_import()
loadfonts(device = "win")
library(ggplot2)
library(readr)
library(dplyr)
data <- read.csv("./Strict/figure1_1.csv", header = FALSE)



names(data)<-c("Cancer","Sex","Percentage")

# Grouped
ggplot(data, aes(fill=Sex, y=Percentage, x=Cancer)) + 
  geom_bar(position="dodge", stat="identity")
the_order<-data$Cancer

png("figure1_1.png", units="cm", width=14, height=8, res=300)
  data %>% 
  ggplot(aes(x = reorder(Cancer, -Percentage), y =Percentage , group = Sex, fill = Sex)) +
  geom_bar(stat = "identity", width = 0.75) +
  coord_flip() +
#  scale_x_discrete(limits = the_order) +
  # another trick!
  scale_y_continuous(breaks = seq(-50, 50, 10), 
                     labels = abs(seq(-50, 50, 10))) +
  labs(x = "Cancer type", y = "Percentage of incidence") +
    theme_classic()+
  theme(legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        text=element_text(size=11, family = "Arial"))+ scale_fill_manual(values=c("#ba0278", "#6aab0f")) 
  dev.off()
  # reverse the order of items in legend
  # guides(fill = guide_legend(reverse = TRUE)) +
  # change the default colors of bars
  #scale_fill_manual(values=c("red", "blue")) 

print(p)
