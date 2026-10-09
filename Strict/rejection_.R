library("readxl")
library(ggplot2)

load("C:/Users/kisl0002/Dropbox/Resistance_2020/Strict/Data/topos_all.RData")
Lapatinib1nm<- read_excel("C:/Users/kisl0002/Downloads/rebound_data/rebound data/Lapatinib1nm_karina_constraint.xls")
Lapatinib10nm<- read_excel("C:/Users/kisl0002/Downloads/rebound_data/rebound data/Lapatinib10nm_karina_constraint.xls")
u0126<- read_excel("C:/Users/kisl0002/Downloads/rebound_data/rebound data/u0126_karina_constraint.xls")
AZD6244<- read_excel("C:/Users/kisl0002/Downloads/rebound_data/rebound data/AZD6244_karina_constraint.xls")
Rapamycin<- read_excel("C:/Users/kisl0002/Downloads/rebound_data/rebound data/Rapamycin_karina_constraint.xls")
Rapamycin<-Rapamycin[-259,]


names(Lapatinib1nm)<-c("Gene","Phosposite","Control","1 min", "5 min", "10 min", "60 min") 
names(Lapatinib10nm)<-c("Gene","Phosposite","Control","1 min", "5 min", "10 min", "60 min") 
names(u0126)<-c("Gene","Phosposite","Control","2.5 min", "5 min", "10 min") 
names(AZD6244)<-c("Gene","Phosposite","Control","15 min", "30 min")
names(Rapamycin)<-c("Gene","Phosposite","Control","2 min", "7 min", "15 min","30 min")


#load("./Strict/Data/network_thesis_withPathways.RData")
load("./Strict/Data/network_thesis.RData")
Lapatinib1nm <- Lapatinib1nm%>%dplyr::distinct(Gene, .keep_all = TRUE)
Lapatinib10nm <- Lapatinib10nm%>%dplyr::distinct(Gene, .keep_all = TRUE)
u0126 <- u0126%>%dplyr::distinct(Gene, .keep_all = TRUE)
AZD6244 <- AZD6244%>%dplyr::distinct(Gene, .keep_all = TRUE)
Rapamycin <- Rapamycin%>%dplyr::distinct(Gene, .keep_all = TRUE)

#################Lapatinib 1nm

nrow(Lapatinib1nm)-length(which(is.element(Lapatinib1nm$Gene,nodes)))
length(which(is.element(Lapatinib1nm$Gene,nodes)))

load("./Strict/Data/results_EGFR_5.RData")
egfr<-results
egfr$Topology <-prodlim::row.match(egfr[,1:9], topos_all)
egfr<-egfr[is.element(egfr$Topology, frec$Topology),] 
egfr_filter<- egfr%>%dplyr::distinct(V12, .keep_all = TRUE)

load("./Strict/Data/results_ERBB2_5.RData")
ebb2<-results
ebb2$Topology <-prodlim::row.match(ebb2[,1:9], topos_all)
ebb2<-ebb2[is.element(ebb2$Topology, frec$Topology),] 
ebb2_filter<- ebb2%>%dplyr::distinct(V12, .keep_all = TRUE)

egfr_ebb2<-unique(c(egfr_filter$V12, ebb2_filter$V12))
inter_lapi1<-Lapatinib1nm[is.element(Lapatinib1nm$Gene,egfr_ebb2),]
length(egfr_ebb2) - nrow(inter_lapi)

Lapatinib1nm_filtered<-Lapatinib1nm[which(is.element(Lapatinib1nm$Gene,egfr_ebb2)),]

#################Lapatinib 10nm

nrow(Lapatinib10nm)-length(which(is.element(Lapatinib10nm$Gene,nodes)))
length(which(is.element(Lapatinib10nm$Gene,nodes)))

inter_lapi10<-Lapatinib10nm[is.element(Lapatinib10nm$Gene,egfr_ebb2),]
length(egfr_ebb2) - nrow(inter_lapi10)

Lapatinib10nm_filtered<-Lapatinib10nm[which(is.element(Lapatinib10nm$Gene,egfr_ebb2)),]

#################AZD6244 10nm

nrow(AZD6244)-length(which(is.element(AZD6244$Gene,nodes)))
length(which(is.element(AZD6244$Gene,nodes)))

load("./Strict/Data/results_MAP2k1_5.RData")
mek1<-results

mek1$Topology <-prodlim::row.match(mek1[,1:9], topos_all)
mek1<-mek1[is.element(mek1$Topology, frec$Topology),] 
mek1_filter<- mek1%>%dplyr::distinct(V12, .keep_all = TRUE)

inter_AZD62<-AZD6244[is.element(AZD6244$Gene,mek1_filter$V12),]
length(mek1_filter$V12) - nrow(inter_AZD62)
nrow(inter_AZD62)

AZD6244_filtered<-AZD6244[which(is.element(AZD6244$Gene,mek1_filter$V12)),]

####################### Rapamycin


nrow(Rapamycin)-length(which(is.element(Rapamycin$Gene,nodes)))
length(which(is.element(Rapamycin$Gene,nodes)))

load("./Strict/Data/results_mTORC1_5.RData")
mtor<-results

mtor$Topology <-prodlim::row.match(mtor[,1:9], topos_all)
mtor<-mtor[is.element(mtor$Topology, frec$Topology),] 
mtor_filter<- mtor%>%dplyr::distinct(V12, .keep_all = TRUE)

inter_Rapamycin<-Rapamycin[is.element(Rapamycin$Gene,mtor_filter$V12),]
length(mtor_filter$V12) - nrow(inter_Rapamycin)
nrow(inter_Rapamycin)





length(which(is.element(Lapatinib1nm$Gene,nodes)))
which(!is.element(Lapatinib1nm$Gene,nodes))




protNames <- read.delim("~/netscan_Beta/data/proteins/uniprot-filtered-organism__Homo+sapiens+(Human)+[9606]_.tab")


found_Lapatinib1nm<-0
found_Lapatinib10nm<-0
found_u0126<-0
found_AZD6244<-0
found_Rapamycin<-0

  for(j in 1:nrow(protNames)){
    genes<-unlist(stringr::str_split(protNames[j,"Gene.names"], pattern =" "))
    
    for(i in 1:nrow(Lapatinib1nm)){
    if(sum(genes== Lapatinib1nm[i,"Gene"])>0){
      found<-c(found_Lapatinib1nm,i)
      break
    }}
      
    for(i in 1:nrow(Lapatinib10nm)){
    if(sum(genes== Lapatinib10nm[i,"Gene"])>0){
      found<-c(found_Lapatinib10nm,i)
      break
    }}
    for(i in 1:nrow(u0126)){
    if(sum(genes== u0126[i,"Gene"])>0){
      found<-c(found_u0126,i)
      break
    }}
    for(i in 1:nrow(AZD6244)){
    if(sum(genes== AZD6244[i,"Gene"])>0){
      found<-c(found_AZD6244,i)
      break
    }}
    for(i in 1:nrow(Rapamycin)){
    if(sum(genes== Rapamycin[i,"Gene"])>0){
      found<-c(found_Rapamycin,i)
      break
    }}
    
    }


to_plot<-AZD6244_filtered

interest_melt<-reshape2::melt(to_plot, id.vars = c("Gene","Phosposite"))
 
#                    measure.vars = c("Control","1 Hour","8 Hour")), color="#005dcf"
#interest_melt<-cbind(id=as.factor(interest_melt[,1]),interest_melt[,-1]), color="#005dcf"

cc <- scales::seq_gradient_pal("#64872c", "#64872c", "Lab")(seq(0,1,length.out=nrow(to_plot)))

png("./Strict/Plots/AZD6244.png", units="cm", width=14, height=7, res=300)
ggplot2::ggplot(interest_melt, aes(x=variable, y= value,color=Gene))+geom_point(size=1)+
  geom_line(aes(group=interaction(Gene,Phosposite), color=Gene), size=0.5) + 
  theme(text=element_text(size=12,colour = "black", family="Heveltica"),
        panel.border = element_blank(),
        #panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        panel.grid = element_line(colour="lightgray", size=0.1), 
        panel.background =  element_rect(fill = "white"),
        axis.text.y= element_text(size=12),
        axis.line.y = element_line(size=1),
        axis.line.x = element_blank(),
        panel.grid.major.x = element_line(colour="black", size=0.5,linetype = 'dashed'),
        legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = element_blank()
        # panel.grid.minor.y = element_line(colour="black", size=1.5),
  )+ ylab("Expression changes") +   scale_colour_manual(values=cc)

#+scale_color_viridis(discrete = TRUE, option = "")


#+ scale_color_manual(values = c("#ebc446", "#41a336", "#FC4E07"))
dev.off()








]
protNames <- read.delim("~/netscan_Beta/data/proteins/uniprot-filtered-organism__Homo+sapiens+(Human)+[9606]_.tab")
load("~/netscan_Beta/data/topologies/topologies_3nodes.RData")
functions_names<- list.files("./data/functions", full.names = T )
lapply(functions_names, source)


interest <- roger[is.element(roger$Proteins, protNames$Entry), ]
interest_2 <- roger[!is.element(roger$Proteins, protNames$Entry), ]


for(j in 1:nrow(interest_2)){
  tmp<- unlist(stringr::str_split(interest_2$Proteins[j], pattern =";"))
  for (k in 1:length(tmp)){
    row_tmp <- interest_2[j,]
    row_tmp$Proteins[1]<-unlist(stringr::str_split(tmp[k], pattern = "-"))[1]
    interest <- rbind(interest,row_tmp)
  }
}

interest_3 <- interest[is.element(interest$Proteins, protNames$Entry), ]


interest_3$gene<-"none"
for (i in 1:nrow(interest_3)){
  interest_3$gene[i] <- protNames[is.element(protNames$Entry, interest_3$Proteins[i]), "Gene.names"]}


roger_rebound<-interest_3[interest_3$cluster=="rebound",]
roger_sensitive<-interest_3[interest_3$cluster=="decrease",]
roger_increase<-interest_3[interest_3$cluster=="increase",]
roger_adaptation<-interest_3[interest_3$cluster=="adaptation",]

to_plot<-roger_rebound

interest_melt<-melt(to_plot, id.vars = c("GPX id","Proteins","Protein names",
                                         "Positions within proteins","Cluster Group","cluster","colorG",
                                         "gene","Phospho STY Probabilities"), 
                    measure.vars = c("Control","1 Hour","8 Hour"))
interest_melt<-cbind(id=as.factor(interest_melt[,1]),interest_melt[,-1])

png("./roger_adaptation.png", units="cm", width=14, height=7, res=300)
ggplot2::ggplot(interest_melt, aes(x=variable, y= value,color=colorG))+geom_point(size=1)+
  geom_line(aes(color=colorG, group=id), size=0.5)+ 
  theme(text=element_text(size=12,colour = "black", family="Heveltica"),
        panel.border = element_blank(),
        #panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        panel.grid = element_line(colour="lightgray", size=0.1), 
        panel.background =  element_rect(fill = "white"),
        axis.text.y= element_text(size=12),
        axis.line.y = element_line(size=1),
        axis.line.x = element_blank(),
        panel.grid.major.x = element_line(colour="black", size=0.5,linetype = 'dashed'),
        legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = element_blank()
        # panel.grid.minor.y = element_line(colour="black", size=1.5),
  )+ ylab("Expression changes") + scale_color_manual(values = c("#1f7f79","#520168", "#E7B800", "#0092e7", "#005dcf"))

#+scale_color_viridis(discrete = TRUE, option = "")


#+ scale_color_manual(values = c("#ebc446", "#41a336", "#FC4E07"))
dev.off()

png("./roger_rebound.png", units="cm", width=14, height=7, res=300)
ggplot2::ggplot(interest_melt, aes(x=variable, y= value))+geom_point(size=1, color="#005dcf")+
  geom_line(aes(group=id), size=0.5, color="#005dcf")+ 
  theme(text=element_text(size=12,colour = "black", family="Heveltica"),
        panel.border = element_blank(),
        #panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        panel.grid = element_line(colour="lightgray", size=0.1), 
        panel.background =  element_rect(fill = "white"),
        axis.text.y= element_text(size=12),
        axis.line.y = element_line(size=1),
        axis.line.x = element_blank(),
        panel.grid.major.x = element_line(colour="black", size=0.5,linetype = 'dashed'),
        legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = element_blank()
        # panel.grid.minor.y = element_line(colour="black", size=1.5),
  )+ ylab("Expression changes") 

#+scale_color_viridis(discrete = TRUE, option = "")


#+ scale_color_manual(values = c("#ebc446", "#41a336", "#FC4E07"))
dev.off()
