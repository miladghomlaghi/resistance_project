
functions_names<- list.files("./Strict/Data/pathways", full.names = T )
pathways<-lapply(functions_names, read.delim)

pathways<-lapply(pathways, function(x){
  x <-x[x$EFFECT != "unknown",]
  x <-x[x$EFFECT != "form complex",]
  x <-x[!is.na(x$EFFECT),]
  x <-x[!is.na(x$SCORE),]
  x <-x[x$DIRECT=="YES",]
  x$EFFECT_BIN<-ifelse(grepl(c("down"),x$EFFECT),-1,1)
  x<-x%>% dplyr::distinct(ENTITYA, ENTITYB, .keep_all = TRUE)
  
  x<- x[,c("ENTITYA","ENTITYB","EFFECT_BIN","EFFECT",
                               colnames(x)[!is.element(colnames(x),
                                                                 c("ENTITYA","ENTITYB","EFFECT_BIN", "EFFECT"))])]
  
 
  x<-x[, c("ENTITYA","ENTITYB","EFFECT_BIN")]
  
  names(x)<-c("ENTITYA","ENTITYB","STATUS")
  return(x)
  })



