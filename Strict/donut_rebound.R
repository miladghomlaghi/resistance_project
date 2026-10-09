library(readxl)

load("./Strict/Data/frequency_matlab_210220.RData")


data <- data.frame(
  category=c("Resistant", "Only functional"),
  count=c(sum(frecuency_Rep$Mean>0),sum(frecuency_Rep$Mean==0))
)



# Compute percentages
data$fraction <- data$count / sum(data$count)

# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)

# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

# Compute a good label
data$label <- paste0(data$count,"\n ",data$category )
cols<-{c("#9aa2b5","#8b008b")}
# Make the plot
png("./Strict/Plots/Circular_proportion.png", units="cm", width=7, height=7, res=300)
ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  #geom_text( x=1.4, aes(y=labelPosition, label=label, color=category), size=12) + # x here controls label position (inner / outer)
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  coord_polar(theta="y") +
  xlim(c(0.5, 4)) +
  theme_void()+
  theme(legend.position = "none")
dev.off()        



