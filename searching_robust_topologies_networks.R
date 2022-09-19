library(foreach)
library(pracma)
library(igraph)
library(doParallel)

load("./Data/Data/minimal_top_thesis.RData")
load("./Data/Data/topos_all.RData")
load("./Data/Data/frequency_matlab_210220.RData")
load("./Data/Data/network_thesis.Rdata")

topos_to_match<- cbind(topos_all[frecb[order(frecb$Mean_BB,decreasing = T),]$Topology[1:8],], rep(7,8),
                       frecb[order(frecb$Mean_BB,decreasing = T),]$Topology[1:8])
names(topos_to_match) <- names(final_good_top)
topos_to_match<- rbind(final_good_top, topos_to_match)


dist<-5
start_node<-5
end_node<-5

pracma::tic()

num_nodes <- 3
positions <- 1:num_nodes^2
p <- 1:num_nodes
nodes<-sort(nodes)
combi_nodes <- gtools::permutations(num_nodes, 2, 1:num_nodes, repeats.allowed=T)
combi_nodes <- combi_nodes[ ,c(2,1)]
range <- 1:length(nodes)

final_results<- data.frame(t(rep("none", 12)))


browser()

for (v in start_node:end_node){
  
  #NL_positions <- num_nodes*(p-1)+p
  
  #user_node_num = which(targets==protein)
  user_node_num <- v
  #positions <- positions[!is.element(positions,NL_positions)]
  
  near_net<-igraph::ego(net,dist, nodes = targets[user_node_num], mode = "out", mindist = 0)
  near_net<-induced_subgraph(net,unlist(near_net))
  
  
  near_nodes<-as_ids(V(near_net))
  
  if(length(near_nodes)>2){
    
    near_name<-which(near_nodes==targets[user_node_num])
    
    all_combi<- gtools::permutations(length(near_nodes)-1, num_nodes-1, 
                                     range[-near_name], repeats.allowed=F)
    #numCores<-parallel::detectCores(all.tests = FALSE, logical = TRUE)
    
    #p1 <- p[-user_position][1]
    #p2 <- p[-user_position][2]
    #god_com <- data.frame(t(rep("none", 11)))
    numCores <- bigstatsr::nb_cores()
    doParallel::registerDoParallel(numCores)
    results <- foreach::foreach(r=1:nrow(all_combi),.combine = "rbind") %dopar% {
      #results <- foreach::foreach(r=5:10,.combine = "rbind") %dopar% {
      #source("./Strict/F05_get_sign.R")
      #source("./Strict/F06_get_link.R")
      source("./R/Strict/F05_get_sign.R")
      source("./R/Strict/F06_get_link.R")
      
      #results <-matrix(rep("none",num_nodes),1,num_nodes)
      
      #for (r in 1:nrow(all_combi)){
      
      combi<- rep(0,num_nodes)
      combi[1] <- near_name
      combi[2]<-all_combi[r,][1]
      combi[3]<-all_combi[r,][2]
      combi <- matrix(combi,nrow=1,ncol = num_nodes)
      combi_names<-near_nodes[combi]############################################
      topo <-rep(0,num_nodes^2)
      
      for (i in 1:nrow(combi_nodes)){
        tmp_nodes <- p[-combi_nodes[i,]] #
        tmp_names <- combi_names[tmp_nodes]
        tmp_net <- igraph::delete_vertices(near_net, tmp_names)
        tmp_link <- get_link(tmp_net,combi_names[combi_nodes[i,1]],combi_names[combi_nodes[i,2]],links)
        if(!is.na(tmp_link)){
          topo[positions[i]]<-tmp_link
        }
      }
      
      # if (!is.na(prodlim::row.match(topo,topos_to_match[,1:9]))){
      #   god_com<-rbind(god_com,c(topo,combi_names))
      # }
      if (sum(topo!=0)>2){
        c(topo, combi_names)
      }
    }
  }
  
  
  doParallel::stopImplicitCluster()
  
  if(exists("results")){
    results<-as.data.frame(results)
    names(results)<-names(final_results)
    final_results<-rbind(final_results,results)
  }
}



#save("final_results", file="./Data/results_Test.RData")

pracma::toc()