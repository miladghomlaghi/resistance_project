find_bad_susp_initial_search <- function(potential_proteins_table,
                                         interaction_matrix,
                                         network_proteins,
                                         all_proteins_permut_proteins,
                                         all_proteins_permut_values,
                                         shorter_col_proteins,
                                         shorter,
                                         longer) {
  
  # browser()
  #####################################################################
  # extracting the rows from the dataset that we will use in this section
  # and remove the rest
  #####################################################################
  
  if (interaction_matrix[shorter, 1] != 0) {
    tmp_index <- intersect(
      which(all_proteins_permut_proteins[, 1] %in% shorter_col_proteins) ,
      which(unlist(all_proteins_permut_proteins[, 2]) == unlist(potential_proteins_table[1, 1]))
    )
    # browser()
    tmp_path_of_proteins_ws <-
      all_proteins_permut_proteins[tmp_index, ]
    
    tmp_path_of_values_ws <-
      all_proteins_permut_values[tmp_index, ]
  }
  
  if (interaction_matrix[1, shorter] != 0) {
    tmp_index2 <- intersect(
      which(all_proteins_permut_proteins[, 2] %in% shorter_col_proteins) ,
        which(unlist(all_proteins_permut_proteins[, 1]) == unlist(potential_proteins_table[1, 1]))
    )
    # browser()
    tmp_path_of_proteins_sw <-
      all_proteins_permut_proteins [tmp_index2,]
    
    
    tmp_path_of_values_sw <-
      all_proteins_permut_values[tmp_index2, ]
  }
  
  rm(all_proteins_permut_proteins, all_proteins_permut_values)
  gc()
  
  ###############################################################################
  # Having tmp_path_value and protein we find suspicious and bad rows for each
  # protein in the shorter protein column
  # #############################################################################
  # browser()
  
  all_susp_bad_rows <-
  foreach::foreach(kk = 1:length(shorter_col_proteins),
                   .combine = "rbind") %dopar% {
                       bad_rows <- NULL
                       next_suspecious_rows <- NULL
                       
                       # for (kk in 1:length(shorter_col_proteins)) { #
                       # if (shorter_col_proteins[kk]==1){
                       #   browser()
                       # }
                       
                       # selecting proteins in the longer column
                       include_starting_index <-
                         which(shorter_col_proteins[kk] == potential_proteins_table[, shorter])
                       longer_col_proteins <-
                         network_proteins[unlist(potential_proteins_table[include_starting_index, longer])]
                       
                       
                       
                       
                       
                       ########################
                       ########################
                       
                       if (interaction_matrix[shorter, 1] != 0) {
                         # finding the paths and values between node 1 and shorter
                         
                         index_tmp_path=which(tmp_path_of_values_ws[,1]==shorter_col_proteins[kk])
                         
                         # rows that dont have protein in the shortest path and their path value is also
                         # different

                           if (interaction_matrix[shorter, 1][[1]] != tmp_path_of_values_ws[[index_tmp_path,3]][1]) {
                           bad_rows <- c(bad_rows, include_starting_index[1] - 1 +
                                           which(
                                             !(longer_col_proteins %in% tmp_path_of_proteins_ws[[index_tmp_path,3]][[1]])
                                           ))
                         }
                         
                         
                         # Finding proteins that are in the 10 paths from node A and B
                         suspecious_rows <-
                           which(longer_col_proteins %in% tmp_path_of_proteins_ws[[index_tmp_path,3]][[1]])
                         
                         
                         if (length(suspecious_rows) != 0) {
                           # merging the list of 10 paths for this interaction
                           path_protein_count <-
                             as.data.frame(table(unlist(tmp_path_of_proteins_ws[[index_tmp_path,3]])))
                           #calculating the accurance rate of the proteins in the 10 paths
                           sample_count_ten <-
                             path_protein_count[path_protein_count[, 1]  %in% longer_col_proteins[suspecious_rows] , 2]
                           
                           # browser()
                           # if a protein exists 10 times in the 10 paths it is not a good candidate
                           if (length(sample_count_ten) != 0)
                             if (any(sample_count_ten == 10)) {
                               bad_rows <-
                                 c(bad_rows,
                                   include_starting_index[1] - 1 + suspecious_rows[sample_count_ten == 10])
                             }
                           next_suspecious_rows <-
                             c(next_suspecious_rows,
                               include_starting_index[1] - 1 + suspecious_rows[sample_count_ten < 10])
                         }
                       }
                       
                       ###########################
                       ###########################
                       
                       if (interaction_matrix[1, shorter] != 0) {
                         # first step is to check if the proteins in longer column are in the
                         # first list of tmp_path_of_proteins_ws
                         
                         # browser()
                         
                         # rows that dont have protein in the shortest path and their path value is also
                         # different
                         index_tmp_path=which(tmp_path_of_values_sw[,2]==shorter_col_proteins[kk])
                         
                         
                         if (interaction_matrix[1, shorter][[1]] != tmp_path_of_values_sw[[index_tmp_path,3]][[1]]) {
                           bad_rows <-
                             c(bad_rows, include_starting_index[1] - 1 + which(
                               !(longer_col_proteins %in% tmp_path_of_proteins_sw[[index_tmp_path,3]][[1]])
                             ))
                         }
                         
                         
                         # Finding proteins that are in the 10 paths from node A and B
                         
                         
                         suspecious_rows <-
                           which(longer_col_proteins %in% tmp_path_of_proteins_sw[[index_tmp_path,3]][[1]])
                         
                         if (length(suspecious_rows) != 0) {
                           path_protein_count <-
                             as.data.frame(table(unlist(tmp_path_of_proteins_sw[[index_tmp_path,3]])))
                           
                           sample_count_ten <-
                             path_protein_count[path_protein_count[, 1] %in%  longer_col_proteins[suspecious_rows] , 2]
                           
                           
                           if (length(sample_count_ten) != 0) {
                             if (any(sample_count_ten == 10)) {
                               bad_rows <-
                                 c(bad_rows,
                                   include_starting_index[1] - 1 + suspecious_rows[sample_count_ten == 10])
                             }
                           }
                           next_suspecious_rows <-
                             c(next_suspecious_rows,
                               include_starting_index[1] - 1 + suspecious_rows[sample_count_ten < 10])
                         }
                         
                       }
                       
                       # all_susp_bad_rows <-
                       list(list(unique(bad_rows)),
                            list(unique(next_suspecious_rows)))
                       
                     }
  
  return(all_susp_bad_rows)
  
}
