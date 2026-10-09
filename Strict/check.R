library('prodlim')
library(dplyr)
library(igraph)

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

load("./Strict/Data/network_thesis.RData")
load("./Strict/Data/minimal_top_thesis.RData")
load("./Strict/Data/frequency_matlab_210220.RData")
load("./Strict/Data/topos_all.RData")

merged_data0$Topology<-prodlim::row.match(merged_data0[,1:9], topos_all)


minimal<-merged_data0[which(is.element(merged_data0$Topology,final_good_top$track)),]

high_robust_top<-merged_data0[which(is.element(merged_data0$Topology,
                                               which(frecuency_Rep$Mean*100>3))),]

high_robust_top<-high_robust_top[which(is.element(high_robust_top$V10,targets)),]
minimal<-minimal[which(is.element(minimal$V10,targets)),]

#high_robust_top_all<-high_robust_top


merged_data_all<-rbind(minimal,high_robust_top)
merged_data<- merged_data_all%>%dplyr::distinct(V10,V12,Topology,.keep_all=T)


targets_freq_all<-table(merged_data_all$V10)
barplot(sort(targets_freq_all))

targets_freq_robust<-as.data.frame(sort(table(high_robust_top$V10),  decreasing = T))
targets_freq_robust$class<-"robust"
targets_freq_minimal<-as.data.frame(sort(table(minimal$V10),  decreasing = T))
targets_freq_minimal$class<-"minimal"

data_bar<- rbind(targets_freq_robust, targets_freq_minimal)

table(data_bar$Var1)


# Small multiple,"#1db3e0","#ffaec8",, "#af38a7""#e01d1d""#3886af","#38af8d"
png("./Strict/Plots/targets_freq.png", units="cm", width=18, height=9, res=300)
ggplot(data_bar, aes(x=Var1, y=Freq, fill=class)) + 
  geom_bar(position="stack", stat="identity") +
 scale_fill_manual(values = c("#38af8d","#af38a7"))+
  # scale_fill_viridis(discrete = T) +
  #ggtitle("Studying 4 species..") +
  theme_classic() +
  #geom_errorbar( aes(x=as.character(Topology), ymin=mean-sd, ymax=mean+sd), width=0.4, colour="gray26", alpha=0.9, size=0.6)+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=10),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10),
    legend.position='none'
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 20000), breaks = seq(0,20000,4000))



dev.off()

high_robust_top_filter<- high_robust_top%>%dplyr::distinct(V10,V12,Topology,.keep_all=T)
minimal_filter<- minimal%>%dplyr::distinct(V10,V12,Topology,.keep_all=T)

targets_freq_robust_filter<-as.data.frame(sort(table(high_robust_top_filter$V10),  decreasing = T))
targets_freq_robust_filter$class<-"robust"
targets_freq_minimal_filter<-as.data.frame(sort(table(minimal_filter$V10),  decreasing = T))
targets_freq_minimal_filter$class<-"minimal"


data_bar<- rbind(targets_freq_robust_filter, targets_freq_minimal_filter)



# Small multiple,"#1db3e0","#ffaec8",, "#af38a7""#e01d1d""#3886af","#38af8d"
png("./Strict/Plots/targets_freq_filter.png", units="cm", width=18, height=9, res=300)
ggplot(data_bar, aes(x=Var1, y=Freq, fill=class)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c("#38af8d","#af38a7"))+
  # scale_fill_viridis(discrete = T) +
  #ggtitle("Studying 4 species..") +
  theme_classic() +
  #geom_errorbar( aes(x=as.character(Topology), ymin=mean-sd, ymax=mean+sd), width=0.4, colour="gray26", alpha=0.9, size=0.6)+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=10),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10),
    legend.position='none'
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1000), breaks = seq(0,1000,100))


dev.off()
unique(data_bar$Var1)



sort(unique(minimal$V12))

outputs<-lapply(high_targets,  function(x,set){return(set[set$V10==x,])}, set=high_robust_top)

high_targets <-unique(high_robust_top_filter$V10)
target_output<-apply(as.matrix(high_robust_top_filter), MARGIN=1, 
                     FUN=function(x,set){return(nrow(set[set$V10==as.character(x[10]) & set$V12==as.character(x[12]) & set$Topology==as.numeric(x[13]),]))}, set=high_robust_top)
target_output<-target_output/max(target_output)
#apply(as.matrix(high_robust_top[1:2,]), MARGIN=1, FUN=function(x,set){print(set$V13)}, set=high_robust_top_all)

#names(outputs)<-high_targets

#ouptpus_freq<-table(high_robust_top$V12)


#barplot(sort(ouptpus_freq[as.logical(ouptpus_freq>2)]))



ouptpus_freq[[6]]
length(unique(high_robust_top$V12))
unique(high_robust_top$V10)

vert<-unique(c(unique(high_robust_top_filter$V12),unique(high_robust_top_filter$V10)))
graph1<-graph_from_data_frame(high_robust_top_filter[,c(10,12)],vertices =vert,directed = F)

graph2 <- visNetwork::toVisNetworkData(graph1)

graph2$nodes$group<-ifelse(is.element(graph2$nodes$id,high_targets), "N","B")
graph2$edges$value <- target_output
c_scale <- colorRamp(c("grey",'darkblue'))
colors_ed = apply(c_scale(target_output), 1, FUN=function(x) {rgb(x[1]/255,x[2]/255,x[3]/255)} )
graph2$edges$color<-colors_ed 

library(visNetwork)
visNetwork(nodes = graph2$nodes ,
           edges = graph2$edges,
           background = 'white', main = NULL)%>%
  visGroups(groupname = "N", color ="red", size=50) %>%
  visGroups(groupname = "B", color ="#3281ea", size=5) %>%
  #visEdges(color = colors_ed)%>%
  #visEdges(arrows = NULL, arrowStrikethrough = F,color = "#3281ea")%>%
  visPhysics(solver = "forceAtlas2Based")


cancer_Genes <- read.delim("./Sequential/Data/2020-02-02_IntOGen-Drivers-20200213/Compendium_Cancer_Genes.tsv", header=T)
cancer_Genes<-unique(cancer_Genes$SYMBOL)

length(unique(high_robust_top[high_robust_top$V10=="SRC",]$V12))
outputs <- high_robust_top[is.element(high_robust_top$V12,cancer_Genes),]
outputs_SRC<-outputs[outputs$V10=="SRC",]
outputs_SRC<-sort(table(outputs_SRC$V12))
outputs_filter <- high_robust_top_filter[is.element(high_robust_top_filter$V12,cancer_Genes),]
ouptpus_freq<-sort(table(outputs$V12))
ouptpus_freq_filter<-sort(table(outputs_filter$V12))

outputs <- minimal[is.element(minimal$V12,cancer_Genes),]
outputs_SMO<-outputs[outputs$V10=="SMO",]
outputs_SMO<-sort(table(outputs_SMO$V12))

length(unique(minimal[minimal$V10=="SMO",]$V12))

outputs <- minimal[is.element(minimal$V12,cancer_Genes),]
outputs_MEK1<-outputs[outputs$V10=="MAP2K1",]
outputs_MEK1<-sort(table(outputs_MEK1$V12))

length(unique(minimal[minimal$V10=="MAP2K1",]$V12))

       