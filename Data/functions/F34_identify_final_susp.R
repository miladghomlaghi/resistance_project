
identify_final_susp <- function(potential_proteins_table,
                                   network_proteins,
                                   all_proteins_permut_proteins,
                                    suspicious_rows,
                                index1){

  all_next_suspicious_rows <- NULL
  all_bad_rows1 <- NULL
  all_bad_rows2 <- NULL
  # browser()
  

if (length(suspicious_rows) != 0) {
  # merging all 10 list of paths for this interaction
  
  
  all_next_suspicious_rows <-
    foreach::foreach(i = 1:length(suspicious_rows),
                     .combine = "rbind") %dopar% {
     
                       
    bad_rows <- NULL
    next_suspicious_rows <- NULL
    
  # for (i in 1:length(suspicious_rows)) {
    
    path_protein_count <-
      as.data.frame(table(unlist(all_proteins_permut_proteins[[i]])))
    sample_count_ten <-
      path_protein_count[path_protein_count[, 1] == network_proteins[potential_proteins_table[1, 1][[1]]] , 2]
    
    if (length(sample_count_ten) != 0){
      if (sample_count_ten == 10) {
        bad_rows <- c(bad_rows,suspicious_rows[[i]])
      }
      
      if (sample_count_ten < 10) {
        next_suspicious_rows <- c(next_suspicious_rows,suspicious_rows[[i]])
        
      }
    }
    list(bad_rows,next_suspicious_rows)
  }
}
  # browser()
  # return(c(list(next_suspicious_rows),list(bad_rows)))
  return(all_next_suspicious_rows)
  # doParallel::stopImplicitCluster()
  gc()
  
}