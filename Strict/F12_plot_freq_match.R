
library(dplyr)

functions_names<- list.files("./Strict/Data/Run_network_thesis/", full.names = T)
load(functions_names[1])
run1<-matching_results[-1,]
load(functions_names[2])
run2<-matching_results[-1,]
load(functions_names[3])
run3<-matching_results[-1,]
load(functions_names[4])
run4<-matching_results[-1,]


merged_data0<-rbind(run1,run2,run3,run4)
merged_data0<-as.data.frame(merged_data0)
load("./Strict/Data/minimal_top_thesis.RData")
load("./Strict/Data/frequency_matlab_210220.RData")
load("./Strict/Data/topos_all.RData")

merged_data0$Topology<-prodlim::row.match(merged_data0[,1:9], topos_all)
#merged_data0<-merged_data0[-1,]

merged_data1<-merged_data0[which(is.element(merged_data0$Topology,final_good_top$track)),]

high_robust_top<-merged_data0[which(is.element(merged_data0$Topology,
                                               which(frecuency_Rep$Mean*100>3))),]


high_robust_top_filter<-high_robust_top%>%dplyr::distinct(V10,V12,Topology,.keep_all=T)
length(unique(high_robust_top_filter$V12))

output_nodes_tab_filter<-table(high_robust_top_filter$V12)
output_input_tab_filter<-table(high_robust_top_filter$V10)
output_input_tab_filter
output_nodes_tab_filter[output_nodes_tab_filter>7]
output_nodes_STAT3 <- high_robust_top_filter[high_robust_top_filter$V10=="SRC","V12"]
unique(output_nodes_STAT3)

table(output_nodes_STAT3)
high_robust_top[is.element(high_robust_top$V12,c("STAT3", "STAT1/STAT3")),]



merged_data<-rbind(merged_data1,high_robust_top)
merged_data<-dplyr::distinct(merged_data)
#############for filtered data
merged_data<- merged_data%>%dplyr::distinct(V10,V12,Topology,.keep_all=T)

minimal_freq<-table(as.factor(merged_data$Topology))


frequency<-as.data.frame(cbind(Individual=as.numeric(names(minimal_freq)),
                 value= ((as.numeric(minimal_freq)*100)/max(as.numeric(minimal_freq)))
                   ))
frequency$group<-"none"

frequency$group[frequency$Individual%in%
                    c(4622,6665,6863,7460 )] <- "NFBL+PFBL"
frequency$group[frequency$Individual%in%
                    c(7820,7841,8213)] <- "IFFL"
frequency$group[frequency$Individual%in%
                  c(8181)] <- "IFFL+NFBL"
frequency$group[frequency$Individual%in%
                  c(8228,10059,10581,15333)] <- "IFFL+PFBL"
frequency$group[frequency$Individual%in%
                    c(6428, 6434, 7925,8009,8039,8258, 9620, 13199, 7616,6260)]<-"NFBL"
frequency$group[frequency$Individual%in%
                  high_robust_top$Topology]<-"High-robust"


table(merged_data$V10)


#frequency_full<-frequency
#frequency <-frequency_full[frequency_full$value>1,]
#frequency_small<-frequency_full[frequency_full$value<1,]



data <- frequency
data <- data %>% dplyr::arrange(group, value)
data$group<-as.factor(data$group)


# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 2
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% dplyr::arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  dplyr::group_by(group) %>% 
  dplyr::summarize(start=min(id), end=max(id) - empty_bar) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

png("./Strict/Plots/ciruclarFreqsignor_nounique.png", units="cm", width=20, height=20, res=300)
# Make the plot

ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +  
geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.6) +
  geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +


  ylim(-100,120) +
  theme_minimal() +
  
  theme(
    legend.position = "none", 
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")
 ) +    
  coord_polar() + 
   geom_text(data=label_data, aes(x=id, y=-35, label=Individual, hjust=hjust, color=group), 
             fontface="bold",alpha=0.9, size=4, angle= label_data$angle, inherit.aes = FALSE ) +
   geom_text(data=data, aes(x=id, y=ifelse(value<10,value+15, value+6), label=value*max(minimal_freq)/100, 
                            hjust=label_data$hjust),color="black", fontface="bold",
             alpha=0.9, size=4, angle= label_data$angle, inherit.aes = FALSE ) +
#
#   # Add base line information
   geom_segment(data=base_data, aes(x = start-0.5, y = -5, 
                                    xend = end+0.5, yend = -5, color=group),  alpha=0.6, size=1.5 , 
                inherit.aes = FALSE )

dev.off()



  # geom_text(data=base_data, aes(x = title, y =-30,label=group),hjust=c(1,1,0,0), colour = "black", alpha=0.6, size=2,angle= c(20,45,-60,-20), fontface="bold", inherit.aes = FALSE)
