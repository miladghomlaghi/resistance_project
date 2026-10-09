

# numCores <- bigstatsr::nb_cores()
# doParallel::registerDoParallel(cl <- makeCluster(numCores))


load("protein_interactions_20_length_27_10_2022.RData")

  #
  # importing results of rebound analysis of experimental data
  rebound_data <-
    read.csv(file = "bounceback_src_shared_two.csv",
             header = TRUE)
  
    rebound_list <- rebound_data
    second_col_proteins <- unique(rebound_list[, 2])
    
    third_col_proteins <- unique(rebound_list[, 3])
    node_index_one <- which(nodes == rebound_list[1, 1])
    
    one_index = which(all_proteins_permut_proteins [, 2] ==
                        node_index_one)
    
    
    ###########################################################################
    # Investigating distance between node one and two
    
    node_index_two<-NULL
    node_index_three<-NULL
    
    for (j in 1:length(second_col_proteins)) {
      node_index_two[j] <- which(nodes == second_col_proteins[j])
      
    }
    for (j in 1:length(third_col_proteins)) {
      node_index_three[j] <- which(nodes == third_col_proteins[j])
      
    }
    
    
    tmp_sec_col_proteins_path <-
      all_proteins_permut_values[intersect(one_index,
                                             which(all_proteins_permut_proteins[, 1] %in% node_index_two)), ]
    
    tmp_sec_col_proteins_first_path <-
      sapply(tmp_sec_col_proteins_path, "[[", 1)


    tmp_thrd_col_proteins_path <-
      all_proteins_permut_values[intersect(one_index,
                                             which(all_proteins_permut_proteins[, 1] %in% node_index_three)), ]
    
    tmp_thrd_col_proteins_first_path <-
      sapply(tmp_thrd_col_proteins_path, "[[", 1)

    
    
    
    