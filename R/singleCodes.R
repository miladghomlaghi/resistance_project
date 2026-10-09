#Sinlge codes


### from 3 replicates, format the frequency data frame and save it as RData
frequency_Rep = cbind(frequency59040_R1Seq, frequency59040_R2Seq[,2],frequency59040_R3Seq[,2])
frequency_Rep <- cbind(frequency_Rep,ceiling(rowMeans(frequency_Rep[,2:4])))

names(frequency_Rep) <- c("Topology", "Rep_1", "Rep_2","Rep_3", "Mean")

save(frequency_Rep,file = 'frequency_Seq_Rep.RData')

frecb <-frequency_Rep[frequency_Rep$Mean>0,]
frec <-frecb