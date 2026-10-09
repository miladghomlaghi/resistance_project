merge_matching_data <- function (folder_name, topos_to_match, topos_all, 
                                 net, nodes, links,nodesAreTargets, targets) {

  
 # folder_name <- "4Nodes/Signor_PathwaysPlus"
  
matching_Datafiles <- list.files(path=paste0("Data/matching/",folder_name,collapse =''))
num_nodes <-sqrt(ncol(topos_all))
p <- 1:num_nodes

targets <- nodes[is.element(nodes,targets)]
targets_num<-match(targets, nodes)


links$ENTITYA<-gsub("/",";",links$ENTITYA) 
links$ENTITYB<-gsub("/",";",links$ENTITYB)
links$ENTITYA<-gsub(",",":",links$ENTITYA) 
links$ENTITYB<-gsub(",",":",links$ENTITYB)
nodes<-gsub("/",";",nodes)
nodes<-gsub(",",":",nodes)
targets <- gsub(",",":",targets)
vertex_attr(net)$name<-gsub("/",";",vertex_attr(net)$name)
vertex_attr(net)$name<-gsub(",",":",vertex_attr(net)$name)

topos_ID <- prodlim::row.match(topos_to_match, topos_all)
frequency_match <- data.frame(topos_ID)
frequency_match[,2]<-0
frequency_match[,3]<-''


for(f in 1:length(matching_Datafiles)) {

  #matching_scores<-foreach (f = 1:length(matching_Datafiles),.combine = "cbind") %dopar% {
 # print(matching_Datafiles[f])
  load(paste0("Data/matching/",folder_name,'/',matching_Datafiles[f],collapse =''))
  
  for (V in 1:nrow(final_results)){
  
    if (final_results[V,2] != "none"){
      tmp_freq <- data.frame("tmpCol"=rep(0,nrow(frequency_match)))
      tmp_freq[,2] <- rep('',nrow(frequency_match))
      optionsV1 <- unlist(strsplit(final_results[V,2], ","))
      range <- 1:length(nodes)
      vv <- as.numeric(final_results[V,1])
      
      all_combi<- gtools::permutations(length(nodes)-1, num_nodes-1, range[-vv], repeats.allowed=F) 
      if (length(nodesAreTargets)>1){
        all_combi<-all_combi[is.element(all_combi[,1],targets_num),]
      }
      
      
      numCores<- parallel::detectCores(all.tests = FALSE, logical = TRUE)
      
      doParallel::registerDoParallel(numCores)
      merge_scores<- foreach (C = 1:length(optionsV1),.combine = "cbind")%dopar%{
        
        for (C in 1:length(optionsV1) ){
      
        source("./R/Functions/F02_get_sign.R")
        source("./R/Functions/F03_get_link.R")
        source("./R/Functions/F10_get_simplified_topology.R")
        
        tmp_scores<- data.frame("tmpCol" = rep(0,nrow(frequency_match)))
        tmp_scores[,2] <-rep('',nrow(frequency_match))
        
        tmp<-as.numeric(optionsV1[C])
        
        
        #combi<-c(vv,all_combi[tmp,])
        combi<- replicate(num_nodes,0)
        if (length(nodesAreTargets)>1){
          combi[nodesAreTargets[1]] <- vv
          combi[nodesAreTargets[2]] <- all_combi[tmp,1]
          combi[p[!is.element(p,nodesAreTargets)]] <- all_combi[tmp,2:(num_nodes-1)] 
        }else{
          combi[nodesAreTargets[1]] <- vv
          combi[p[!is.element(p,nodesAreTargets)]] <- all_combi[tmp,3:(num_nodes-1)]
        }
        
        combi_names <- nodes[combi]
        
        
        topo <- get_simplified_topology(combi_names, net, links)
        
        tmp_scores[prodlim::row.match(topo,topos_to_match),1]<-
          tmp_scores[prodlim::row.match(topo,topos_to_match),1]+1
        
        tmp_scores[prodlim::row.match(topo,topos_to_match),2]<-
          paste0(tmp_scores[prodlim::row.match(topo,topos_to_match),2],
                 paste0(combi_names, collapse = ','),collapse = '/')
        
         # print(tmp_scores[prodlim::row.match(topo,topos_to_match),2],V)
        
      #  data.frame(tmp_scores)
     # }
      
      
     
    #   stopImplicitCluster()
    #   
    #   if(ncol(merge_scores)!=2){
    #     tmp_freq[,1] <- rowSums(merge_scores[,seq(1,ncol(merge_scores),2)])
    #     tmp_freq[,2] <- apply(merge_scores[,seq(2,ncol(merge_scores),2)],1, 
    #                           function(row) paste(row[nzchar(row)], collapse = "/")
    #     )}else{
    #       tmp_freq[,1] <- merge_scores[,1]
    #       tmp_freq[,2] <- merge_scores[,2]
    #     }
    #   frequency_match[,2] <- rowSums(cbind(frequency_match[,2],tmp_freq[,1]))
    #   frequency_match[,3] <- apply(cbind(frequency_match[,3],tmp_freq[,2]),
    #                                1,function(row) paste(row[nzchar(row)], collapse = "/"))
    # }
    
}
}

#return(frequency_match)
        return(topo)
}






