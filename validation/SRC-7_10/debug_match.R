## debug_match.R  (self-contained; as.character guarantees no stray column names)
## Run AFTER net, nodes, links, targets, protein_index, topos_to_match,
## combi_nodes, positions, num_nodes, p exist.

library(dplyr)
source("../../R/Strict/F05_get_sign.R")
source("../../R/Strict/F06_get_link.R")
source("../../R/Strict/F07_get_link_path.R")

match_one <- function(r, all_combi, near_name, near_nodes, near_net,
                      combi_nodes, positions, topos_to_match,
                      num_nodes, p, links) {
  
  combi <- rep(0, num_nodes)
  combi[1] <- near_name
  combi[2] <- all_combi[r, ][[1]]
  combi[3] <- all_combi[r, ][[2]]
  combi <- matrix(combi, nrow = 1, ncol = num_nodes)
  ## as.character() returns a PLAIN unnamed character vector.
  cn <- as.character(near_nodes[unlist(combi)])
  
  topo  <- rep(0, num_nodes ^ 2)
  paths <- vector("list", num_nodes ^ 2)
  
  for (i in 1:nrow(combi_nodes)) {
    tmp_nodes <- p[-combi_nodes[i, ]]
    tmp_names <- cn[tmp_nodes]
    tmp_net   <- igraph::delete_vertices(near_net, tmp_names)
    tmp_res   <- get_link_path(tmp_net, cn[combi_nodes[i, 1]],
                               cn[combi_nodes[i, 2]], links)
    if (!is.na(tmp_res$sign)) {
      topo[positions[i]]    <- tmp_res$sign
      paths[[positions[i]]] <- tmp_res$path
    }
  }
  
  row_list <- list()
  for (ii in 1:nrow(topos_to_match)) {
    nz <- as.vector(topos_to_match[ii, ] != 0)
    if (!is.na(prodlim::row.match(topo[nz],
                                  as.data.frame(topos_to_match[ii, nz]), nomatch = NA))) {
      for (pos in which(nz)) {
        fr  <- combi_nodes[pos, 1]
        to  <- combi_nodes[pos, 2]
        pth <- paths[[pos]]
        row_list[[length(row_list) + 1]] <- data.frame(
          X1 = cn[1], X2 = cn[2], X3 = cn[3],
          topo_id = as.integer(ii),
          highly_robust = ii >= 35,
          edge_pos = as.integer(pos),
          from_slot = as.integer(fr), to_slot = as.integer(to),
          from_name = cn[fr], to_name = cn[to],
          sign = as.numeric(topo[pos]),
          path = if (length(pth) >= 1) paste(pth, collapse = ">") else "",
          path_len = if (length(pth) >= 1) as.integer(length(pth) - 1) else NA_integer_,
          n_intermediate = as.integer(max(0, length(pth) - 2)),
          intermediates = if (length(pth) > 2)
            paste(pth[2:(length(pth) - 1)], collapse = ">") else "",
          stringsAsFactors = FALSE)
      }
    }
  }
  
  if (length(row_list) == 0) return(NULL)
  dplyr::bind_rows(row_list)
}

stopifnot(grepl("as.character\\(near_nodes", paste(deparse(body(match_one)), collapse = " ")))

near_net   <- net
near_nodes <- nodes
range      <- 1:length(near_nodes)
near_name  <- which(near_nodes == "SRC")

all_combi <- expand.grid(range[-near_name], protein_index)
all_combi <- all_combi[-which(all_combi$Var1 == all_combi$Var2), ]
cat("Total combinations:", nrow(all_combi), "\n")

expected_cols <- c("X1","X2","X3","topo_id","highly_robust","edge_pos",
                   "from_slot","to_slot","from_name","to_name","sign","path",
                   "path_len","n_intermediate","intermediates")

N <- min(300, nrow(all_combi))
bad <- c(); collected <- list()
for (r in 1:N) {
  out <- tryCatch(
    match_one(r, all_combi, near_name, near_nodes, near_net,
              combi_nodes, positions, topos_to_match, num_nodes, p, links),
    error = function(e) { cat(sprintf("  [r=%d] ERROR: %s\n", r, conditionMessage(e)))
      structure("ERROR", class = "match_error") })
  if (inherits(out, "match_error")) { bad <- c(bad, r); next }
  if (is.null(out)) next
  if (!identical(names(out), expected_cols)) {
    cat(sprintf("  [r=%d] COLUMN MISMATCH:\n", r)); print(names(out)); bad <- c(bad, r)
  }
  collected[[length(collected) + 1]] <- out
}

cat("\nmatched chunks:", length(collected), " | problem rows:", length(bad), "\n")
if (!length(bad)) cat("All matched chunks have the correct 15 columns.\n")

if (length(collected)) {
  combined <- dplyr::bind_rows(collected)
  cat("bind_rows OK. rows:", nrow(combined), " cols:", ncol(combined), "(want 15)\n")
  print(utils::head(combined))
}