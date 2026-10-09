#### Merge pathway from signor

```{r}
signor_EGFR <- read.delim("./RowData/Databases/Signor/Pathways/signor_EGFR.tsv")
signor_Hippo <- read.delim("./RowData/Databases/Signor/Pathways/signor_Hippo_signalling.tsv")
signor_insulin <- read.delim("./RowData/Databases/Signor/Pathways/signor_insulin.tsv")

signor_mtor <- read.delim("./RowData/Databases/Signor/Pathways/signor_mtor.tsv")
signor_p38 <- read.delim("./RowData/Databases/Signor/Pathways/signor_p38.tsv")
signor_PIK3AKT <- read.delim("./RowData/Databases/Signor/Pathways/signor_PIK3AKT.tsv")


signor_pathways <- rbind(signor_EGFR,signor_Hippo,signor_insulin,signor_mtor,signor_p38,signor_PIK3AKT)

signor_pathways$ENTITYA<-as.character(signor_pathways$ENTITYA)
signor_pathways$ENTITYB<-as.character(signor_pathways$ENTITYB)



signor_pathways <-signor_pathways[signor_pathways$EFFECT != "unknown",]
signor_pathways <-signor_pathways[signor_pathways$EFFECT != "form complex",]

#signor_pathways$EFFECT<-ifelse(grepl(c("down"),signor_pathways[,'EFFECT']),-1,1)
signor_pathways <-signor_pathways[signor_pathways$DIRECT != "NO",]

signor_pathways <- signor_pathways[signor_pathways$TAX_ID==9606,]

#signor_pathways <- signor_pathways%>%distinct(ENTITYA, ENTITYB,EFFECT, .keep_all = TRUE)
####signor_pathways <- signor_pathways%>%distinct(ENTITYA, ENTITYB, PMID, .keep_all = TRUE)

signor_pathways<-signor_pathways[complete.cases(signor_pathways[,c("ENTITYA","ENTITYB")]),]

nodesSignorPathways<-unique(c(signor_pathways$ENTITYA,signor_pathways$ENTITYB))
linksSignorPathways<-signor_pathways[,-which(names(signor_pathways)%in% c("ENTITYA","ENTITYB","EFFECT"))]
linksSignorPathways<-cbind(signor_pathways[,which(names(signor_pathways)%in% c("ENTITYA","ENTITYB","EFFECT"))],linksSignorPathways)
netSignorPathways= graph_from_data_frame(linksSignorPathways,vertices = nodesSignorPathways,directed = T)
save(linksSignorPathways,nodesSignorPathways,netSignorPathways,file = 'Signor_Pathways.Rdata')

signor_pathways_plus<-linksSignorPathways
#a=2870

for (a in 1:nrow(signor_filter)){
  if (signor_filter$ENTITYA[a]%in%nodesSignorPathways & signor_filter$ENTITYB[a]%in%nodesSignorPathways  ){
    if(is.na( row.match(c(signor_filter$ENTITYA[a],signor_filter$ENTITYB[a]),linksSignorPathways[,c(1,2)]))){
      signor_pathways_plus<- rbind(signor_pathways_plus,signor_filter[a,] )}
  }
}
signor_pathways_plus <- signor_pathways_plus%>%distinct(ENTITYA, ENTITYB,EFFECT, .keep_all = TRUE)
signor_pathways_plus<-signor_pathways_plus[complete.cases(signor_pathways_plus[,c("ENTITYA","ENTITYB")]),]




nodesSignorPathwaysPlus<-unique(c(signor_pathways_plus$ENTITYA,signor_pathways_plus$ENTITYB))
linksSignorPathwaysPlus<-signor_pathways_plus[,-which(names(signor_pathways_plus)%in% c("ENTITYA","ENTITYB","EFFECT"))]
linksSignorPathwaysPlus<-cbind(signor_pathways_plus[,which(names(signor_pathways_plus)%in% c("ENTITYA","ENTITYB","EFFECT"))],linksSignorPathwaysPlus)
netSignorPathwaysPlus= graph_from_data_frame(linksSignorPathwaysPlus,vertices = nodesSignorPathwaysPlus,directed = T)
save(linksSignorPathwaysPlus,nodesSignorPathwaysPlus,netSignorPathwaysPlus,file = 'Signor_PathwaysPlus.Rdata')


```