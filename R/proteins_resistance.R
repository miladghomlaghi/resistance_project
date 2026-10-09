
library(stringr)

protein_action <- read.delim("C:/Users/kisl1/Downloads/Signalling_Files/9606.protein.actions.v11.0.txt")
protein.aliases <- read.delim("R_Data/String/9606.protein.aliases.v11.0.txt", header=FALSE, comment.char="#")
BIOGRID <- read.delim("C:/Users/kisl1/Downloads/PhD/Resistance-CodeKarina/R_codes/Biogrid/BIOGRID-ALL-3.5.181.tab.txt",header=T, sep="\t")
ACSN2 <- read.delim("C:/Users/kisl1/Downloads/Signalling_Files/ACSN2/ACSN2_binary_relations_between_proteins_with_PMID.txt")
drugs_FDA <- read_excel("C:/Users/kisl1/Downloads/Signalling_Files/drugs_FDA.xlsx", 
                        +     col_names = FALSE)

names(drugs_FDA)<- c("Drug","Area","Target")
names(protein.aliases)<- c("NodeString","Name","Origin")

protein.aliases <- protein.aliases[grep("BLAST_KEGG_NAME",protein.aliases$Origin),]

library("methods")
data <- read_xml('C:/Users/kisl1/Downloads/Signalling_Files/KEEG/hsa05200_Cancer.XML')
entry <- xml_find_all(data, ".//entry")
relation <- xml_find_all(data, ".//relation")

relation_df<-data.frame(matrix(,length(relation), ncol=5))
for(i in 1:length(relation)){
  tmp1<-stringr::str_extract(string = relation[i], pattern = "entry1=\"[0-9]+")
  tmp2<-stringr::str_extract(string = relation[i], pattern = "entry2=\"[0-9]+")
  tmpTy<-stringr::str_extract(string = relation[i], pattern = "type=\"[A-z]+")
  tmpSt<-stringr::str_extract(string = relation[i], pattern = "subtype name=\"[A-z]+")
  tmpVa<-stringr::str_extract(string = relation[i], pattern = "value=\".+\"")
  tmp1<-substr(tmp1,9,nchar(tmp1))
  tmp2<-substr(tmp2,9,nchar(tmp2))
  tmpTy<-substr(tmpTy,7,nchar(tmpTy))
  tmpSt<-substr(tmpSt,15,nchar(tmpSt))
  tmpVa<-substr(tmpVa,8,nchar(tmpVa)-1)
  relation_df[i,]<-cbind(tmp1,tmp2,tmpTy,tmpSt,tmpVa)[1,]
}

names(relation_df)<-c("entry1","entry2","TypeRelation","TypeRelationName","TypeRelationSymbol")

entry_df<-data.frame(matrix(,length(entry), ncol=5))
for(i in 1:length(entry)){
  tmpID<-stringr::str_extract(string = entry[i], pattern = "entry id=\"[0-9]+")
  tmpNa<-stringr::str_extract(string = entry[i], pattern = "name=\".+\" t")
  tmpTy<-stringr::str_extract(string = entry[i], pattern = "type=\"[A-z]+")
  tmpGr<-stringr::str_extract(string = entry[i], pattern = "graphics name=\".+\" f")
  tmpCo<-stringr::str_extract_all(string = entry[i], pattern = "component id=\"[0-9]+")
  tmpCo<-paste(stringr::str_extract_all(string = tmpCo, pattern = "[0-9]+")[[1]], collapse = ",")
  tmpID<-substr(tmpID,11,nchar(tmpID))
  tmpNa<-substr(tmpNa,7,nchar(tmpNa)-3)
  tmpTy<-substr(tmpTy,7,nchar(tmpTy))
  tmpGr<-substr(tmpGr,16,nchar(tmpGr)-3)
  entry_df[i,]<-cbind(tmpID,tmpNa,tmpTy,tmpGr,tmpCo)[1,]
}

names(entry_df)<-c("ID","Code","TypeElement","NameElemnt","Component")
entry_df$NameElemnt<-str_remove(entry_df$NameElemnt,"\\.+")
entry_df$NameElemnt<-str_remove(entry_df$NameElemnt,"\\+/-+")
entry_df$NameElemnt<-str_remove(entry_df$NameElemnt,"\\+/-+")
entry_df$NameElemnt<-str_remove(entry_df$NameElemnt,"\\(")
entry_df$NameElemnt<-str_remove(entry_df$NameElemnt,"\\)")

for (i in 332:nrow(entry_df)){
  tmpList <-str_split(entry_df[i,5], ",")
  finalName<-""
  for(j in 1:length(tmpList[[1]])){
    entry_df[tmpList[[1]][j],4]
    finalName<-paste(finalName,entry_df[tmpList[[1]][j],4],sep ="/")
  }
  entry_df[i,4]<-paste("GR:",substr(finalName,2, nchar(finalName)))
}

entry_df$New_ID<-""

unique_element<-unique(entry_df$NameElemnt)[-1]

for (i in 1:length(unique_element)){
  entry_df$New_ID[grep(unique_element[i], entry_df$NameElemnt)]<-paste(grep(unique_element[i], entry_df$NameElemnt)[1],"a",sep="")
}



entry_relation<-relation_df
for(i in 1:nrow(entry_relation)){
  entry_relation$entry1[i]<-entry_df$New_ID[entry_relation$entry1[i]==entry_df$ID ]
  entry_relation$entry2[i]<-entry_df$New_ID[entry_relation$entry2[i]==entry_df$ID ]
}

#grep(str_split(drugs_FDA$Target[2],"/")[[1]][1],entry_relation$entry1)

entry_df$Name_FDA<-""
drugs_FDA$Target_Entry<-""


for(i in 1:nrow(drugs_FDA)){
  flag = 1
  target<- str_trim(str_split(drugs_FDA$Target[i],"/")[[1]])
  for (j in 1:length(target)){
    if(flag>0){
     matches<-grep(target[j],entry_df$NameElemnt)
    if(length(matches)>0) { 
    for (k in 1:length(matches)){
      matchList<- str_remove(str_trim(str_split(entry_df$NameElemnt[matches[k]], ",")[[1]]),"\\.+")
      #match(FALSE,matchList==target)
      if(sum(matchList==target[j])>0){
        entry_df$Name_FDA[matches[k]]<-target[j]
        drugs_FDA$Target_Entry[i]<-target[j]
        flag = 0
      }}}}}}


subBioGrid <- BIOGRID[,c(8,9,20)]

subBioGrid <- subBioGrid[is.element(subBioGrid$Modification,c("Phosphorylation","Dephosphorylation")),] 


for(i in 1:nrow(BIOGRID)){
  inter1<- str_trim(BIOGRID[i,8])
  inter2<- str_trim(BIOGRID[i,9])
  match1<-grep(inter1,entry_df$NameElemnt)
  match2<-grep(inter2,entry_df$NameElemnt)
  if(length(match1)>0 & length(match2)>0){
    tmp1 <- 0
    for(j in 1:length(match1)){
      if (sum(inter1==str_trim(str_split(entry_df$NameElemnt[match1[j]], ",")[[1]]))==1){
        tmp1<-1}}
    tmp2<-0
    for(j in 1:length(match2)){
      if (sum(inter2==str_trim(str_split(entry_df$NameElemnt[match2[j]], ",")[[1]]))==1){
        tmp2<-1}}
    if (tmp1+tmp2>1){
    subBioGrid<- rbind(subBioGrid,BIOGRID[i,c(8,9,12,13,20)])}
  }
}
subBioGrid<-subBioGrid[-1,]
 
subBioGrid$Inter1<-""
subBioGrid$Inter2<-""
for(i in 1:nrow(subBioGrid)){
  inter1<- str_trim(subBioGrid[i,1])
  inter2<- str_trim(subBioGrid[i,2])
  match1<-grep(inter1,entry_df$NameElemnt)
  match2<-grep(inter2,entry_df$NameElemnt)
  subBioGrid$Inter1[i]<- entry_df$New_ID[match1[1]]
  subBioGrid$Inter2[i]<- entry_df$New_ID[match2[1]]
}
