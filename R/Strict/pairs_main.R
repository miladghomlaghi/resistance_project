#### Plot correlation between parameters

library(readxl)
library(ggplot2)
library(GGally)
library(hrbrthemes)
library(fields)
library(reshape2)
library(extrafont)
library(dplyr)
library(ggpubr)

top_1="7841"
top_2="8214"

kcats <- read_excel(paste0("Strict/Data/Ks/kcats_",top_2,"_",top_1,".xls",collapse=''),col_names = FALSE)
kms<- read_excel(paste0("Strict/Data/Ks/kms_",top_2,"_",top_1,".xls",collapse=''),col_names = FALSE)
names_ks <- read_excel(paste0("Strict/Data/Ks/names_",top_2,"_",top_1,".xls",collapse=''),col_names = FALSE)
sets_ks <- read_excel(paste0("Strict/Data/Ks/sets_",top_2,"_",top_1,".xls",collapse=''),col_names = FALSE)
sets_ks_main <- read_excel(paste0("Strict/Data/Ks/sets_",top_1,".xls",collapse=''),col_names = FALSE)
sets_ks_extra <- read_excel(paste0("Strict/Data/Ks/sets_",top_2,".xls",collapse=''),col_names = FALSE)

compare<-rep("both",length(sets_ks))
compare[is.element(sets_ks,sets_ks_main) & !is.element(sets_ks,sets_ks_extra)]<-paste("T",top_1,sep='')
compare[is.element(sets_ks,sets_ks_extra) & !is.element(sets_ks,sets_ks_main)]<-paste("T",top_2,sep='')
compare
common_sets<-sets_ks[compare=="both"]
sum(compare==paste("T",top_2,sep=''))
sum(compare==paste("T",top_1,sep=''))

ks<- log10(kcats/kms)/3
names(ks)<-names_ks[1,]
ks$compare<-compare
ks$compare<-as.factor(ks$compare)













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
    geom_point(alpha=0.9)+
    scale_x_continuous(limits = c(-1,1),
                       breaks = c(-1,-0.5,0,0.5,1), labels = c("-1","","0","","1"))+
    scale_y_continuous(limits = c(-1,1),
                       breaks = c(-1,-0.5,0,0.5,1), labels = c("-1","","0","","1"))
    #scale_color_viridis(discrete = TRUE, alpha=0.7)     
}  

three_col<-c("#c7028f","#fbb41a","#400967")

three_col<-c("#f28e13","#fbb41a","#400967")
png(paste("Strict/Plots/pairs_main_",top_2,"_",top_1,".png",sep = ''), units="cm", width=24, height=16, res=300)
ggpairs(ks,columns = c(1:7),mapping = ggplot2::aes(colour=compare),
       
      upper = list(continuous = wrap(upperfun)),
      # upper = list(continuous = wrap(ggally_cor, size=3.5, digits=2)),
       # upper = list(continuous = wrap(ggally_cor, size=3.5, digits=)),
        lower = list(continuous = wrap(lowerfun)),
        diag = list(continuous = wrap("densityDiag", alpha=0.8))) +theme_bw() +
  # scale_color_manual(values = c("#7AAE3D", "#7554A3", "#f58a02"))+
  # scale_fill_manual(values = c("#7AAE3D", "#7554A3", "#f58a02"))
  scale_color_manual(values =three_col)+
  scale_fill_manual(values = three_col) +scale_x_continuous(limits = c(-1,1),
                                                           breaks = c(-1,-0.5,0,0.5,1), labels = c("-1","","0","","1"))
# scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
# scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
dev.off()


top_1="7841"
kcats <- read_excel(paste0("Strict/Data/Ks/kcats_",top_1,".xls",collapse=''),col_names = FALSE)
kms<- read_excel(paste0("Strict/Data/Ks/kms_",top_1,".xls",collapse=''),col_names = FALSE)
names_ks <- read_excel(paste0("Strict/Data/Ks/names_",top_1,".xls",collapse=''),col_names = FALSE)
sets_ks<- read_excel(paste0("Strict/Data/Ks/sets_",top_1,".xls",collapse=''),col_names = FALSE)

overshoot<- c(206,	898,	1335,	1341,	2142,	2259,	2872,	2920,	3134,
             3538,	3674,	3747,	3803,	3900,	4100,	4142,	4165,	4239,	4323,	4373,	4760,	4813,	4990)

overshoot<-c(1341,	1513,	2371,	2513,	3323)
 

ks<- log10(kcats/kms)/3
names(ks)<-names_ks[1,]
ks$compare<-paste0("T",top_1,collapse='')
ks$compare<-as.factor(ks$compare)

ks$over<-"rebound"
ks[match(overshoot,sets_ks),"over"]<-"over"
ks$over<-as.factor(ks$over)

sum(ks$compare=="both")
sum(ks$compare==paste("T", top_1,sep=''))


color_dots<- "#80c9b3"
color_dots<-"#da641a"
png(paste("Strict/Plots/pairs_over_",top_1,".png",sep = ''), units="cm", width=14, height=10, res=300)
ggpairs(ks,columns = c(2,3,4,5,6,7),mapping = ggplot2::aes(colour=over),
        
        #upper = list(continuous = wrap(upperfun)),
        upper = list(continuous = wrap(ggally_cor, size=3.5, digits=2)),
        # upper = list(continuous = wrap(ggally_cor, size=3.5, digits=)),
        lower = list(continuous = wrap(lowerfun)),
        diag = list(continuous = wrap("densityDiag", alpha=0.8))) +theme_bw() +
  # scale_color_manual(values = c("#7AAE3D", "#7554A3", "#f58a02"))+
  # scale_fill_manual(values = c("#7AAE3D", "#7554A3", "#f58a02"))+
  scale_color_manual(values =three_col)+
  scale_fill_manual(values =three_col)+
  scale_x_continuous(limits = c(-1,1), breaks = c(-1,-0.5,0,0.5,1), labels = c("-1","","0","","1"))
# scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
# scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
dev.off()



png(paste("Strict/Plots/pairs_main_",top_1,".png",sep = ''), units="cm", width=24, height=16, res=300)
ggpairs(ks_2,columns = c(2,4,5,6),mapping = ggplot2::aes(colour=compare),
        
        upper = list(continuous = wrap(upperfun)),
        #upper = list(continuous = wrap(ggally_cor, size=3.5, digits=2)),
        # upper = list(continuous = wrap(ggally_cor, size=3.5, digits=)),
        lower = list(continuous = wrap(lowerfun)),
        diag = list(continuous = wrap("densityDiag", alpha=0.8))) +theme_bw() +
  # scale_color_manual(values = c("#7AAE3D", "#7554A3", "#f58a02"))+
  # scale_fill_manual(values = c("#7AAE3D", "#7554A3", "#f58a02"))
  scale_color_manual(values =three_col)+
  scale_fill_manual(values = three_col) +scale_x_continuous(limits = c(-1,1),
                                                            breaks = c(-1,-0.5,0,0.5,1), labels = c("-1","","0","","1"))
# scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
# scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
dev.off()
