## F07_get_link_path.R
## New function: like get_link() but also returns the node path (intermediates).
## The SIGN logic is identical to get_link(), so topology matching is unchanged.
## Returns a list: $sign (numeric, as get_link) and $path (character vector of
## node names from protA to protB, or NA when no path exists).

get_link_path <- function(network, protA, protB, links_frame) {

  if (protA != protB) {

    path <- names(unlist(igraph::shortest_paths(network, protA, protB)[[1]]))

    if (is.null(path) || length(path) == 0) {
      path_value <- 0
      return(list(sign = path_value, path = NA_character_))
    } else {
      path_value <- 1
      for (p in 1:(length(path) - 1)) {
        tmp <- links_frame[prodlim::row.match(c(path[p], path[p + 1]),
                             as.data.frame(links_frame[, c(1, 2)])), 3]
        path_value <- c(path_value, tmp)
      }
      path_value <- prod(path_value)
      return(list(sign = path_value, path = path))
    }

  } else {
    ## self-loop: same node to same node
    path_value <- links_frame[prodlim::row.match(c(protA, protB),
                               as.data.frame(links_frame[, c(1, 2)])), 3]
    if (length(path_value) == 0) path_value <- NA
    return(list(sign = path_value, path = protA))
  }
}
