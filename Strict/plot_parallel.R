library(viridis)
library(ggplot2)
library("RColorBrewer")
library(GGally)
library("readxl")
top_1 <- "6431"
top_2 <- "8213"
top_3 <- "6434"

path_top<-"Strict/Data/Ks/sets_"


load("./Strict/Data/frequency_matlab_210220.RData")
data_all<-frecuency_Rep
data_all[as.numeric(top_1),]
data_all[as.numeric(top_2),]
data_all[as.numeric(top_5),]

para_names<-data.frame(t(c(rep(0,9),"I","I","e11","e11","e21","e21","e31","e31",
                           "e12","e12","e22","e22","e32","e32","e13","e13",
                           "e23","e23","e33","e33", "B1","B1","B2","B2","B3","B3")))


categrories<-c(paste0(top_1,"-",top_2,"-", top_3),
               paste0(top_1,"-", top_2),
               paste0(top_1,"-", top_3),
               paste0(top_2,"-", top_3),
               paste0("Only", top_1),
               paste0("Only", top_2),
               paste0("Only", top_3))
groups<-1:7
               

sets_T1<- read_excel(paste0(path_top,top_1,".xls", collapse = ""),col_names = FALSE)
sets_T2 <- read_excel(paste0(path_top,top_2,".xls", collapse = ""),col_names = FALSE)
sets_T3 <- read_excel(paste0(path_top,top_3,".xls", collapse = ""),col_names = FALSE)

sets_T1 <- data.frame(set=as.numeric(t(sets_T1)))
sets_T2 <- data.frame(set=as.numeric(t(sets_T2)))
sets_T3 <- data.frame(set=as.numeric(t(sets_T3)))

tmp_sets<-rbind(sets_T1,sets_T2,sets_T3)

set_all<-unique(tmp_sets$set)
paste0(set_all,collapse = ",")
set_all<-data.frame(set=set_all)

set_all$Group<-"none"

for (i in 1:nrow(set_all)){
  if(is.element(set_all[i,1],t(sets_T1))){
    if(is.element(set_all[i,1],t(sets_T2))){
      if(is.element(set_all[i,1],t(sets_T3))){
        set_all[i,2]<-1
      }else{set_all[i,2]<-2}
    }else{
      if(is.element(set_all[i,1],t(sets_T3))){
        set_all[i,2]<-3
      }else{
        set_all[i,2]<-5
      }
    }
  }else{
    if(is.element(set_all[i,1],t(sets_T2))){
      if(is.element(set_all[i,1],sets_T3)){
        set_all[i,2]<-4
      }else{set_all[i,2]<-6}
    }else{
      set_all[i,2]<-7}
    }
  }

set_all$Group<-as.factor(set_all$Group)

set_all$set <- as.factor(set_all$set)
summary(set_all)

set_all <- set_all[order(set_all$set),]


mypalette<-brewer.pal(7,"Accent")
mypalette <- c("#ffd70e","#ff7f27" ,"#ff5fc4","black", "#a00390","#0ed145", "#5f8cff","#696653", "#16bbc9")


kcats_3T <- read_excel(paste0("Strict/Data/Ks/kcats_",top_1,"_",top_2,"_",top_3,".xls",collapse=''),col_names = FALSE)
kms_3T <- read_excel(paste0("Strict/Data/Ks/kms_",top_1,"_",top_2,"_",top_3,".xls",collapse=''),col_names = FALSE)
names_ks <- read_excel(paste0("Strict/Data/Ks/names_",top_1,"_",top_2,"_",top_3,".xls",collapse=''),col_names = FALSE)
sets_ks <- read_excel(paste0("Strict/Data/Ks/sets_",top_1,"_",top_2,"_",top_3,".xls",collapse=''),col_names = FALSE)

ks <- log10(kcats_3T/kms_3T)/3
#ks_s<-log10(kcats_s/kms_s)/3
#ks_s<-ks_s[-3]
names(ks) <- names_ks

ks$ID <- as.factor(t(sets_ks))

ks <- ks[order(ks$ID),]
#names(ks_s)<-names_ks[1,]
#ks_s$compare<-"special"
ks$compare <- set_all$Group

ks$compare <- as.factor(ks$compare)


a <- 1
png("./Strict/Plots/parallel_fullInter.png", units="cm", width=14, height=7, res=300)
ggparcoord(ks[is.element(ks$compare,a),], scale = "globalminmax",
           columns = c(1:8), groupColumn = 10, 
           showPoints = TRUE, 
           title = "Parallel Coordinate Plot for the Iris Data",
           alphaLines =1
          # mapping=aes(size=0.5)
) +
  #scale_color_viridis(discrete=TRUE) +
  scale_color_manual(values=c(mypalette[a],mypalette[a]))+
  #theme_bw()+
  theme(text=element_text(size=12,colour = "black", family="Arial"),
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
  ) + scale_y_continuous(limits = c(-1,1), sec.axis = dup_axis()) + 
  ylab("Regulation strength") 

dev.off()


set_all$set[set_all$Group==6]

ks<-ks[!ks$ID==1704,]
a<-3
b<-7 
c<-6
ks[is.element(ks$compare,c(a,b)) & ks$Stimulus<0,"ID"]
ks[ks$compare==7 & ks$Input>0,"ID"]
png("./Strict/Plots/parallel_green_orange_less0.png", units="cm", width=14, height=7, res=300)
ggparcoord(ks[is.element(ks$compare,c(a,b)) & ks$Input<0,], scale = "globalminmax",
           #columns = c(1:(c-1),(c+1):7), groupColumn = 10, 
           columns = c(1:8), groupColumn = 10,
           showPoints = TRUE, 
           title = "Parallel Coordinate Plot for the Iris Data",
           alphaLines =1
           # mapping=aes(size=0.5)
) +
  #scale_color_viridis(discrete=TRUE) +
  scale_color_manual(values=c(mypalette[a],mypalette[b]))+
  #theme_bw()+
  theme(text=element_text(size=12,colour = "black", family="Arial"),
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
  ) + scale_y_continuous(limits = c(-1,1), sec.axis = dup_axis()) + 
  ylab("Regulation strength") 

dev.off()


ks<-ks[ks$ID==1704,]
a<-3
b<- 7 
ks[ks$compare==6& ks$Stimulus<0,"ID"]
png("./Strict/Plots/parallel_two_extra_sl0.png", units="cm", width=14, height=7, res=300)
ggparcoord(ks[is.element(ks$compare,c(a,b))& ks$Input<0,], scale = "globalminmax",
           columns = c(1:7), groupColumn = 10, 
           showPoints = TRUE, 
           title = "Parallel Coordinate Plot for the Iris Data",
           alphaLines =1
           # mapping=aes(size=0.5)
) +
  #scale_color_viridis(discrete=TRUE) +
  scale_color_manual(values=c(mypalette[a],mypalette[b]))+
  #theme_bw()+
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
  ) + scale_y_continuous(limits = c(-1,1), sec.axis = dup_axis()) + 
  ylab("Kinetic efficiency") 

dev.off()


a=6
s=140
png("./Strict/Plots/parallel_s78.png", units="cm", width=14, height=7, res=300)
ggparcoord(ks[which(ks$ID==s),], scale = "globalminmax",
           columns = c(1:8), groupColumn = 10, 
           showPoints = TRUE, 
           title = "Parallel Coordinate Plot for the Iris Data",
           alphaLines =1,
           
           # mapping=aes(size=0.5)
) +
  scale_color_viridis(discrete=TRUE) +
  scale_color_manual(values=c(mypalette[a],mypalette[a]))+
  #theme_bw()+
  theme(text=element_text(size=12,colour = "black", family="Arial"),
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
  ) + scale_y_continuous(limits = c(-1,1), sec.axis = dup_axis()) + 
  ylab("Regulation strength") 
dev.off()












n<-7
nrow(ks[is.element(ks$compare, n),])
categrories[n]

paste0(ks[is.element(ks$compare, c(categrories[n])),]$ID, collapse = ",")
bad_T8213<-ks[ks$compare==categrories[n] &  is.element(ks$ID, c(3747,3622,3523,502,2259,3460,1911,1301,2089,1341,898,830,1027)),]
good_T8213<-ks[ks$compare==categrories[n] &  !is.element(ks$ID, c(3747,3622,3523,502,2259,3460,1911,1301,2089,1341,898,830,1027)),]

ks[is.element(ks$compare, c(categrories[n]))&ks$BA>0 & ks$AB<0,]


sets_T3$set

upperfun <- function(data,mapping){
  ggplot(data = data, mapping = mapping)
  #   #geom_smooth() +
  #   #    ylim(0,1.1e5) +
  #   stat_cor(aes(
  #                label =paste(..r.label.., cut(..p.., 
  #                                               breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
  #                                               labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), 
  #                             sep = "~")), 
  #          label.x = 1)
  #label.x = -1, label.y = 2) + 
  #+
  # coord_cartesian(ylim = c(0,1.25), xlim = c(0,0.7))
  #geom_density2d()+
  #scale_x_continuous(limits = c(-1,1))+
  #scale_y_continuous(limits = c(-1,1))
}   

lowerfun <- function(data,mapping){
  ggplot(data = data, mapping = mapping)+
    geom_point(alpha=0.6)+
    scale_x_continuous(limits = c(-1,1),
                       breaks = c(-1,-0.5,0,0.5,1), labels = c("-1","","0","","1"))+
    scale_y_continuous(limits = c(-1,1),
                       breaks = c(-1,-0.5,0,0.5,1), labels = c("-1","","0","","1"))
  #scale_color_viridis(discrete = TRUE, alpha=0.7)     
}  


three_col<-c("#400967","#c7028f","#fbb41a", "black","blue","gray")
png(paste0("cor_",top_3,".png", collapse = ""), units="cm", width=24, height=16, res=300)
ggpairs(ks[is.element(ks$ID,sets_T1$set),],columns = 1:8,mapping = ggplot2::aes(colour="#400967"),
          upper = list(continuous = wrap(upperfun)),
        #upper = list(continuous = wrap(ggally_cor, size=3.5, digits=2)),
        # upper = list(continuous = wrap(ggally_cor, size=3.5, digits=)),
        lower = list(continuous = wrap(lowerfun)),
        diag = list(continuous = wrap("densityDiag", alpha=0.8))) +theme_bw() +
  # scale_color_manual(values = c("#7AAE3D", "#7554A3", "#f58a02"))+
  # scale_fill_manual(values = c("#7AAE3D", "#7554A3", "#f58a02"))
  scale_color_manual(values =three_col)+
  scale_fill_manual(values = three_col) +  scale_x_continuous(limits = c(-1,1),
                                                              breaks = c(-1,-0.5,0,0.5,1), labels = c("-1","","0","","1"))
# scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
# scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
dev.off()

