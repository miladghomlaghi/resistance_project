## x2_distance_analysis.R
## For each rebounded X3, how far is the X2 node from BOTH ends:
##   dist_X1X2 = length of the SRC(X1) -> X2 path  (edge from_slot=1, to_slot=2)
##   dist_X2X3 = length of the X2 -> X3 path       (edge from_slot=2, to_slot=3)
## Reduced per X3 by the X2 that is CLOSEST to X3 (min dist_X2X3), reporting
## both arms for that bridge, plus the separate minimum of each arm.
##
## Reads edge_level_SRC.csv (100MB) via data.table::fread.
## Outputs:
##   x3_x2_distance_highlyrobust.csv  -> question 1 (highly robust only)
##   x3_x2_distance_all.csv           -> question 2 (all matched topologies)

library(data.table)

edge <- data.table::fread("edge_level_SRC.csv")

## --- the two arms ---
e12 <- edge[from_slot == 1 & to_slot == 2,
            .(X2, X3, topo_id, highly_robust,
              dist_X1X2 = path_len, path_X1X2 = path)]
e23 <- edge[from_slot == 2 & to_slot == 3,
            .(X2, X3, topo_id, highly_robust,
              dist_X2X3 = path_len, path_X2X3 = path)]

## join the two arms of the same set (same X2, X3, topology)
bridge <- merge(e12, e23, by = c("X2", "X3", "topo_id", "highly_robust"),
                all = FALSE)

## ---------- helper: per-X3 summary over its X2 bridges ----------
summarise_by_x3 <- function(dt) {
  ## collapse to one row per (X2,X3): an X2 may serve several topos; keep the
  ## arm lengths (identical across topos for the same shortest paths)
  b <- unique(dt[, .(X2, X3, dist_X1X2, dist_X2X3, path_X1X2, path_X2X3)])
  b[, .(
    n_X2_bridges   = uniqueN(X2),
    ## closest bridge to X3 (primary reducer = min X2->X3 distance)
    min_X2X3       = min(dist_X2X3),
    closest_X2     = X2[which.min(dist_X2X3)],
    its_X1X2       = dist_X1X2[which.min(dist_X2X3)],   # SRC->X2 arm of that bridge
    path_X2X3      = path_X2X3[which.min(dist_X2X3)],
    path_X1X2      = path_X1X2[which.min(dist_X2X3)],
    ## separate minima of each arm (best case on each side, any bridge)
    min_X1X2       = min(dist_X1X2),
    median_X2X3    = as.numeric(median(dist_X2X3)),
    median_X1X2    = as.numeric(median(dist_X1X2))
  ), by = X3][order(min_X2X3, X3)]
}

## ---------- Question 1: highly robust only ----------
q1 <- summarise_by_x3(bridge[highly_robust == TRUE])
data.table::fwrite(q1, "x3_x2_distance_highlyrobust.csv")

## ---------- Question 2: all matched ----------
q2 <- summarise_by_x3(bridge)
data.table::fwrite(q2, "x3_x2_distance_all.csv")

## ---------- console summary ----------
cat("Q1 (highly robust):", nrow(q1), "rebounded X3 proteins\n")
cat("   closest X2->X3 distance (per X3):\n");  print(table(q1$min_X2X3))
cat("   closest SRC->X2 distance (per X3):\n"); print(table(q1$min_X1X2))

cat("\nQ2 (all matched):", nrow(q2), "rebounded X3 proteins\n")
cat("   closest X2->X3 distance (per X3):\n");  print(table(q2$min_X2X3))
cat("   closest SRC->X2 distance (per X3):\n"); print(table(q2$min_X1X2))

cat("\nHead of Q1:\n"); print(utils::head(q1))
cat("\nHead of Q2:\n"); print(utils::head(q2))