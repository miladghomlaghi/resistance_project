## Libraries used across the scripts
library(gtools)
library(colorspace)
library(dendextend)#
library(visNetwork)#
library(DiagrammeR)#
library(DiagrammeRsvg)#
library(dplyr)#
library(foreach)#
library(ggplot2)#
library(GGally)#
library(ggExtra)#
library(gridExtra)
library(magrittr)
library(network)
library(png)
library(prodlim)#
library(readxl) #
library(writexl)
library(reshape2)#
library(RColorBrewer)#
library(rsvg)#
library(rstudioapi)
library(sna)
library(scales)
library(tidyverse)
library(svglite)#
library(igraph)#
library(doParallel)
library(ggpubr)
library(ggthemes)
library(gridExtra)
####

library(gstat)
library(windfarmGA)
library(viridis)
library(numbers)
####
#library(BBmisc)#

# install.packages('pracma')
# install.packages('colorspace')
# install.packages('dendextend')#
# install.packages('visNetwork')#
# install.packages('DiagrammeR')#
# install.packages('DiagrammeRsvg')#
# install.packages('dplyr')#
# install.packages('foreach')#
# install.packages('ggplot2')#
# install.packages('GGally')#
# install.packages('ggExtra')#
# install.packages('gridExtra')
# install.packages('magrittr')
# install.packages('network')
# install.packages('png')
# install.packages('prodlim')#
# install.packages('readxl') #
# install.packages('"xlsx"')
# install.packages('reshape2')#
# install.packages('RColorBrewer')#
# install.packages('rsvg')#
# install.packages('rstudioapi')
# install.packages('sna')
# install.packages('scales')
# install.packages('tidyverse')
# install.packages('svglite')#
# install.packages('igraph')#
# install.packages('doParallel')
# ####

# install.packages('gstat')
# install.packages('windfarmGA')
# install.packages('viridis')
# install.packages('numbers')
# ####
# install.packages('BBmisc')





#Data to load
#load('R_Data/topologies_Resistance_Data_12Dec.RData')

# #Load the files with data from Matlab
# frecuency_Rep <- read_excel("freq_resistanceFilter_V2_2.xls", col_names =FALSE)
# over_Rep <- read_excel("overShoot_Filter_V2_2.xls", col_names =FALSE)
# perceValuePlus <- read_excel("perceRobust_resistance3SPlus_V2_2.xls",col_names = FALSE)
# perceValue <- read_excel("perceRobust_resistance3S_V2_2.xls", col_names =FALSE)
# topos <- read_excel("topos_3Nodes.xls",col_names = FALSE)
# 
# #Calculate mean and max of replicates
# foreach (j = 1:nrow(frecuency_Rep)) %do%  (frecuency_Rep[j,11]<-mean(as.numeric(frecuency_Rep[j,c(3,6,9)])))
# foreach (j = 1:nrow(frecuency_Rep)) %do%  (frecuency_Rep[j,12]<-max(as.numeric(frecuency_Rep[j,c(3,6,9)])))
# foreach (j = 1:nrow(frecuency_Rep)) %do%  (frecuency_Rep[j,13]<-mean(as.numeric(frecuency_Rep[j,c(4,7,10)])))
# foreach (j = 1:nrow(frecuency_Rep)) %do%  (frecuency_Rep[j,14]<-max(as.numeric(frecuency_Rep[j,c(4,7,10)])))
# foreach (j = 1:nrow(frecuency_Rep)) %do%  (frecuency_Rep[j,15]<-mean(as.numeric(frecuency_Rep[j,c(2,5,8)])))
# foreach (j = 1:nrow(frecuency_Rep)) %do%  (frecuency_Rep[j,16]<-max(as.numeric(frecuency_Rep[j,c(2,5,8)])))
# 
# #Change the names 
# names(frecuency_Rep) <- c("Topology", "Rep_1_BB","Rep_1_GS","Rep_1","Rep_2_BB","Rep_2_GS","Rep_2",
#                           "Rep_3_BB","Rep_3_GS","Rep_3", "Mean_GS","Max_GS","Mean","Max","Mean_BB","Max_BB")
# names(over_Rep) <- c("Topology", "Rep_1_BB","Rep_1_OV","Rep_1","Rep_2_BB","Rep_2_OV","Rep_2",
#                      "Rep_3_BB","Rep_3_OV","Rep_3", "Mean_OV", "Mean_BB","Mean_Ratio")


# #Only the Mean of replicates, removing those topologies where there were not functional sets if any.
# frecGS<-frecuency_Rep[,c(1,11,13,15)]
# frecGS[is.na(frecGS)]<-0
# frecGS<-frecGS[!frecGS$Mean_GS==0,]
# 
# #Remove topologies where there were any bounce back sets
# frec<-frecuency_Rep
# frec[is.na(frec)]<-0
# frec <-frec[!frec$Mean==0,]
# 
# #Subset roboust topologies
# frecb<-frecuency_Rep
# frecb[is.na(frecb)]<-0
# frecb <-frecb[frecb$Mean>0.01,]
# 
# #Remove from overshoot those topologies where there was not bounce back
# over<-over_Rep
# over[is.na(over)]<-0
# over$Mean<-frecuency_Rep$Mean
# over <-over[!over$Mean_BB ==0,]
# 
# 
# #Melt data frames
# meltedFre <- melt(frec[,c("Topology","Mean")],id.vars = 1) #all topologies
# meltedFre$value <- meltedFre$value *100
#meltedFreb <- melt(frecb[c("Topology","Rep_1", "Rep_2","Rep_3","Mean")],id.vars = 1) #topo more than 100, replicates and mean
# meltedFreb$value <- meltedFreb$value *100
# meltedFre2b <- melt(frecb[,c("Topology","Mean")],id.vars = 1) #topo more thatn 100, only mean

#save.image("topologies_Resistance_Data_12Dec.RData")


#install.packages("GGally")
#install.packages("network")
#install.packages("sna")
#install.packages("DiagrammeR")
#install.packages("readxl")
#install.packages("gplots")
#install.packages("ggdendro")
#install.packages("dendextend")
#install.packages("colorspace")
#install.packages("ggExtra")
#install.packages("prodlim")
#install.packages("magrittr")
#install.packages("rsvg")
#install.packages("DiagrammeRsvg")
#install.packages("svglite")
