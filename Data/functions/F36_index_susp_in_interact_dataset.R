index_susp_in_interact_dataset <- function(interaction_matrix,
                                                   potential_proteins_table,
                                                   all_proteins_permut_values,
                                                   all_next_suspicious_rows,
                                                   shorter,
                                                   longer) {
  # finding index of the suspicous rows in all_proteins_permut_values
  index1 <- NULL
  index2 <- NULL
  
  if (interaction_matrix[shorter, longer] != 0) {
    index1 <-
      row.match(unique(potential_proteins_table[all_next_suspicious_rows, c(shorter, longer)]),
                all_proteins_permut_values[, c(1, 2)])
    
  }
  
  if (interaction_matrix[longer, shorter] != 0) {
    index2 <-
      row.match(unique(potential_proteins_table[all_next_suspicious_rows, c(longer, shorter)]),
                all_proteins_permut_values[, c(1, 2)])
  }
  index <- unique(c(index1, index2))

return(index)
}

