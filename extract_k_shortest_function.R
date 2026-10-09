path_proteins_function <-
  function(net,
           nodes,
           links,
           tmp_proteins_permut) {
    path_proteins <- list()
    
    path_proteins <-
    foreach::foreach(r = 1:nrow(tmp_proteins_permut),
                     .combine = "rbind") %dopar% {
                         
                         # for(r in 1:nrow(tmp_proteins_permut)){
                           
                         source("Data/functions/F02_get_link_topology_extract.R")
                         source("Data/functions/F19_calculate_path_value.R")
                         # print(r)
                         path_values_proteins <- list()
                         
                         
                         path_values_proteins <-
                           get_link_topology_extract(net,
                                                     nodes[[tmp_proteins_permut[r, 1]]],
                                                     nodes[[tmp_proteins_permut[r, 2]]],
                                                     links,
                                                     full = F)
                         
                         list(tmp_proteins_permut[r, 1],
                              tmp_proteins_permut[r, 2],
                              path_values_proteins[[1]],
                              path_values_proteins[[2]])
                         
                         
                       }
    gc()

    return(path_proteins)
    
  }