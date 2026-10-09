## src_to_x3_distance.R
## Distance from SRC to each rebounded X3, THROUGH the motif bridge:
##   dist_SRC_X3 = dist(SRC->X2) + dist(X2->X3)   for a given (X2, topology)
## Reduced per X3 to the MINIMUM over all X2 bridges (shortest motif route).
## Reported for (1) highly robust motifs and (2) all matched motifs.

library(data.table)

edge <- data.table::fread("edge_level_SRC.csv")

e12 <- edge[from_slot == 1 & to_slot == 2,
            .(X2, X3, topo_id, highly_robust, d12 = path_len)]
e23 <- edge[from_slot == 2 & to_slot == 3,
            .(X2, X3, topo_id, highly_robust, d23 = path_len)]

## same set = same X2, X3, topology: add the two arms
b <- merge(e12, e23, by = c("X2","X3","topo_id","highly_robust"), all = FALSE)
b[, dist_SRC_X3 := d12 + d23]

## per-X3 shortest motif route, for a given topology subset
dist_table <- function(dt) {
  per_x3 <- dt[, .(min_SRC_X3 = min(dist_SRC_X3)), by = X3]
  list(n = nrow(per_x3), tab = table(per_x3$min_SRC_X3), per_x3 = per_x3)
}

q1 <- dist_table(b[highly_robust == TRUE])   # highly robust
q2 <- dist_table(b)                          # all matched

cat("Highly robust:", q1$n, "rebounded X3 proteins\n")
cat("SRC->X3 distance distribution (shortest motif route):\n")
print(q1$tab)
cat("\nAll matched:", q2$n, "rebounded X3 proteins\n")
cat("SRC->X3 distance distribution (shortest motif route):\n")
print(q2$tab)

## save per-protein distances
data.table::fwrite(q1$per_x3, "src_to_x3_highlyrobust.csv")
data.table::fwrite(q2$per_x3, "src_to_x3_all.csv")
