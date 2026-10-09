bulk_identify_bad_susp <- function(potential_proteins_table,
         interaction_matrix,
         network_proteins,
         all_proteins_permut_proteins,
         all_proteins_permut_values,
         shorter,
         longer,
         index1){




# extracting the first list of each selected row

first_paths <-
  lapply(all_proteins_permut_proteins[index1, 3], `[[`, 1)

# check if in each row has the user protein
# browser()
tmp_3_bad_rows <- which(unlist(lapply(all_proteins_permut_values[index1, 3], `[[`, 1)) != interaction_matrix[shorter, longer][[1]])

suspicious_rows <-
  which(unlist((
    lapply(first_paths, function(x)
      network_proteins[potential_proteins_table[1, 1][[1]]] %in% x)
  )))


second_suspicious_rows <- 1:nrow(potential_proteins_table)
second_suspicious_rows <- second_suspicious_rows[-suspicious_rows]


all_bad_rows1 <- intersect(tmp_3_bad_rows,second_suspicious_rows)
# browser()
return(c(list(suspicious_rows),list(all_bad_rows1)))
}