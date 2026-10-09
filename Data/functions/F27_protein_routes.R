protein_routes <- function(network_proteins, links, network) {
  all_proteins_permut <-
    gtools::permutations(length(network_proteins),
                         2,
                         1:length(network_proteins),
                         repeats.allowed = F)
  # load("permut.rData")
  
  check_link_exist = 0
  path_values <- list()
  path_proteins <- list()
  # browser()
  
  path_proteins <-
    foreach::foreach(r = 1:nrow(all_proteins_permut),
                     .combine = "rbind") %dopar% {
                       source("./data/functions/F02_get_link_topology_extract.R")
                       source("./data/functions/F19_calculate_path_value.R")
                       
                       # for (r in 1:nrow(all_proteins_permut[1:100,])) {
                       print(r)
                       path_values_proteins <-
                         get_link_topology_extract(network,
                                                   network_proteins[all_proteins_permut[r, 1]],
                                                   network_proteins[all_proteins_permut[r, 2]],
                                                   links,
                                                   full = F)
                       
                       list(path_values_proteins[[1]], path_values_proteins[[2]])
                       
                       
                     }
  all_proteins_permut_proteins = cbind(all_proteins_permut, path_proteins[, 2])
  all_proteins_permut_values = cbind(all_proteins_permut,  path_proteins[, 1])
  return(list(all_proteins_permut_proteins, all_proteins_permut_values))
}