## check_predictions.R
## Take NetScan predicted motifs (X1, X2, X3 per row, highly-robust, no length
## limit) and flag which remain valid under a distance rule.
##
## Needs `net` (the directed igraph SIGNOR network) already loaded, same object
## the pipeline uses. Distances are directed, out (SRC acts downstream), hop
## counts = number of edges, exactly matching the `path_len` used earlier.

library(igraph)
library(data.table)

## ---- 1. read the predictions ------------------------------------------------
pred <- data.table::fread("C:/Users/milad/Desktop/resistance_09_12_final/Results/length_7_phospho_bind_7_6/SRC.csv")   # <-- your file name

## adapt if needed: assume first three columns are X1, X2, X3
setnames(pred, 1:3, c("X1", "X2", "X3"))
pred[, X1 := as.character(X1)]
pred[, X2 := as.character(X2)]
pred[, X3 := as.character(X3)]

## ---- 2. distance helper (directed, out) -------------------------------------
## igraph::distances returns the hop count, or Inf if unreachable.
V_names <- igraph::V(net)$name
in_net  <- function(x) x %in% V_names

dist_out <- function(a, b) {
  if (!in_net(a) || !in_net(b)) return(NA_real_)   # node not in network
  d <- igraph::distances(net, v = a, to = b, mode = "out")[1, 1]
  if (is.infinite(d)) NA_real_ else d               # unreachable -> NA
}

## ---- 3. compute the three arm distances per row -----------------------------
## vectorise over unique node pairs to avoid recomputing
uniq_pairs <- unique(rbind(
  pred[, .(a = X1, b = X2)],
  pred[, .(a = X2, b = X3)],
  pred[, .(a = X1, b = X3)]
))

np <- nrow(uniq_pairs)
cat("Computing distances for", np, "unique node pairs",
    "(", nrow(pred), "prediction rows)...\n")
t0  <- Sys.time()
d   <- numeric(np)
pb  <- utils::txtProgressBar(min = 0, max = np, style = 3)
step <- max(1, floor(np / 200))   # update the bar ~200 times
for (i in seq_len(np)) {
  d[i] <- dist_out(uniq_pairs$a[i], uniq_pairs$b[i])
  if (i %% step == 0 || i == np) {
    utils::setTxtProgressBar(pb, i)
    el  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    eta <- el / i * (np - i)
    cat(sprintf("  | %d/%d (%.0f%%)  elapsed %.0fs  eta %.0fs   \r",
                i, np, 100 * i / np, el, eta))
  }
}
close(pb)
cat("\nDistances done in", round(as.numeric(difftime(Sys.time(), t0,
                                                     units = "secs"))), "s\n")

uniq_pairs[, d := d]
setkey(uniq_pairs, a, b)

getd <- function(a, b) uniq_pairs[.(a, b), d]
pred[, d12 := getd(X1, X2)]
pred[, d23 := getd(X2, X3)]
pred[, d13 := getd(X1, X3)]

## ---- 4. apply the rule ------------------------------------------------------
## CHANGE THIS LINE to test any rule:
pred[, valid := (d12 == 1) & (d23 <= 3)]
## e.g. stricter:  pred[, valid := (d12 == 1) & (d23 <= 3) & (d13 <= 4)]
## rows with NA distance (node missing / unreachable) are NOT valid
pred[is.na(valid), valid := FALSE]

## ---- 5. output + summary ----------------------------------------------------
data.table::fwrite(pred, "netscan_predictions_checked.csv")

cat("Total predicted motif rows:        ", nrow(pred), "\n")
cat("Rows still valid under the rule:    ", sum(pred$valid), "\n")
cat("Reduction:                          ",
    round(100 * (1 - sum(pred$valid) / nrow(pred)), 1), "%\n\n")

cat("Unique predicted X3 (node 3) total: ", uniqueN(pred$X3), "\n")
cat("Unique predicted X3 still valid:    ", uniqueN(pred[valid == TRUE]$X3), "\n")

## how many rows dropped because a node was missing/unreachable
cat("\nRows with a missing/unreachable arm:",
    sum(is.na(pred$d12) | is.na(pred$d23)), "\n")