# =============================================================================
# plot_three_node_paths.R
#
# Given THREE proteins (nodes), using the same signed, directed interaction
# network as the topology-search pipeline, this script produces:
#
#   (A) DETAILED path plot  -> three_node_paths.png / .pdf
#       The SHORTEST route between each ordered pair, INCLUDING intermediate
#       ("middle") proteins.
#
#   (B) SIMPLIFIED net-effect topology -> three_node_net_topology.png / .pdf
#       Each chain A -> ... -> B collapsed to ONE edge whose sign is the PRODUCT
#       of the edge signs along the shortest route (= net activation/inhibition).
#       Also prints a 3x3 net-effect matrix (row = from, col = to).
#
# Edge style (both plots):  activation  --->   (filled arrowhead, green)
#                           inhibition  ---|    (flat T-bar,       red)
# Edges are drawn by hand in base R so the ---| inhibition head is correct;
# igraph is used only to lay out / find paths.
#
# Pipeline mapping:
#   F02 get_link_topology_extract(): shortest-first paths  -> we take the shortest route.
#   F19 calculate_path_value(..., full = FALSE): prod() of signs  -> the net effect.
#
# Run from the project root (beside the .RData):
#   Rscript plot_three_node_paths.R MAPK3 AKT3 CHUK
# =============================================================================

# ---------------------------- CONFIG -----------------------------------------
three_nodes  <- c("MAPK3", "AKT3", "CHUK")   # default; override via args
data_file    <- "new_initial_info_9_6_phospho_bind.RData"
out_detailed <- "three_node_paths"
out_simple   <- "three_node_net_topology"
route_mode   <- "shortest"   # "shortest" (single) or "all" (all shortest, detailed plot only)

ACT_COL <- "#2CA25F"; INH_COL <- "#DE2D26"; NA_COL <- "grey55"
MAIN_FILL <- "#E8833A"; INTER_FILL <- "#9ECAE1"
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) three_nodes <- args[1:3]
if (length(three_nodes) != 3) stop("Please supply exactly three protein names.")
if (!requireNamespace("igraph", quietly = TRUE))
  stop("The 'igraph' package is required: install.packages('igraph').")
library(igraph)

# ---- Load network -----------------------------------------------------------
if (!file.exists(data_file)) stop("Cannot find '", data_file, "'. Run from project root.")
load(data_file)                          # nodes, links, targets, net
if (!exists("links")) stop("'links' not found in ", data_file)

edges <- data.frame(from = as.character(links[[1]]),
                    to   = as.character(links[[2]]),
                    sign = suppressWarnings(as.numeric(links[[3]])),
                    stringsAsFactors = FALSE)
edges <- edges[stats::complete.cases(edges[, c("from","to","sign")]), ]
edges <- edges[edges$from != edges$to, ]
g <- igraph::graph_from_data_frame(edges, directed = TRUE)

edge_key <- paste(edges$from, edges$to, sep = "\r")
sign_of  <- function(a, b) edges$sign[match(paste(a, b, sep = "\r"), edge_key)]

present <- three_nodes %in% igraph::V(g)$name
if (!all(present)) stop("Not in the network: ", paste(three_nodes[!present], collapse = ", "))

# ---- Routes between each ordered pair ---------------------------------------
shortest_route <- function(g, a, b) {
  sp <- suppressWarnings(igraph::shortest_paths(g, a, b, mode = "out")$vpath[[1]])
  if (length(sp) < 2) return(NULL)
  igraph::V(g)$name[as.integer(sp)]
}
all_shortest_routes <- function(g, a, b) {
  res <- suppressWarnings(igraph::all_shortest_paths(g, a, b, mode = "out")$res)
  res <- Filter(function(p) length(p) > 1, res)
  lapply(res, function(p) igraph::V(g)$name[as.integer(p)])
}
net_effect <- function(route) {
  s <- sign_of(route[-length(route)], route[-1])
  if (any(is.na(s))) return(NA_real_)
  prod(s)
}

pair_grid <- expand.grid(a = three_nodes, b = three_nodes, stringsAsFactors = FALSE)
pair_grid <- pair_grid[pair_grid$a != pair_grid$b, ]

chosen_routes <- list(); detailed_routes <- list()
net_matrix <- matrix(NA_real_, 3, 3, dimnames = list(three_nodes, three_nodes))
report <- character(0)
for (i in seq_len(nrow(pair_grid))) {
  a <- pair_grid$a[i]; b <- pair_grid$b[i]
  one <- shortest_route(g, a, b)
  if (is.null(one)) { report <- c(report, sprintf("  %s -> %s : no directed path", a, b)); next }
  ne <- net_effect(one); net_matrix[a, b] <- ne
  chosen_routes[[paste(a, b)]] <- one
  mid <- if (length(one) > 2) paste(one[-c(1, length(one))], collapse = " -> ") else "(direct)"
  report <- c(report, sprintf("  %s -> %s : %d hop(s), net = %s | middle: %s",
                              a, b, length(one) - 1,
                              ifelse(is.na(ne), "NA", ifelse(ne > 0, "+ activation", "- inhibition")), mid))
  detailed_routes <- c(detailed_routes,
                       if (route_mode == "all") all_shortest_routes(g, a, b) else list(one))
}
cat("Route summary (shortest route per ordered pair):\n", paste(report, collapse = "\n"),
    "\n\nNet-effect matrix (row=from, col=to; +1 activation, -1 inhibition, NA none):\n", sep = "")
print(net_matrix); cat("\n")
if (length(detailed_routes) == 0) stop("No directed routes between any pair.")

# =============================================================================
# Hand-drawn renderer: activation ---> , inhibition ---|
# =============================================================================
# coords: n x 2 (data units, plotted with asp = 1)
# radii : per-node circle radius (same units)
# e_df  : data.frame(fi, ti, sign)  edges by node index
draw_signed_graph <- function(coords, labels, radii, fills, label_cex, e_df, title,
                              curv, headlen, headwid, barhalf, lwd, pad) {
  xr <- range(coords[,1]); yr <- range(coords[,2])
  plot(NA, xlim = c(xr[1]-pad, xr[2]+pad), ylim = c(yr[1]-pad, yr[2]+pad),
       asp = 1, axes = FALSE, xlab = "", ylab = "", main = title)
  
  heads <- list()  # draw heads after nodes so they sit on top
  # --- pass 1: edge lines (trimmed to node borders) ---
  for (k in seq_len(nrow(e_df))) {
    A <- coords[e_df$fi[k], ]; B <- coords[e_df$ti[k], ]
    rA <- radii[e_df$fi[k]];   rB <- radii[e_df$ti[k]]
    s  <- e_df$sign[k]
    col <- if (is.na(s)) NA_COL else if (s > 0) ACT_COL else INH_COL
    d <- B - A; L <- sqrt(sum(d^2)); if (L == 0) next
    u <- d / L; perp <- c(-u[2], u[1])
    Ctrl <- (A + B)/2 + perp * curv * L               # reciprocal edges curve apart
    tt <- seq(0, 1, length.out = 240)
    x <- (1-tt)^2*A[1] + 2*(1-tt)*tt*Ctrl[1] + tt^2*B[1]
    y <- (1-tt)^2*A[2] + 2*(1-tt)*tt*Ctrl[2] + tt^2*B[2]
    keep <- which(sqrt((x-A[1])^2+(y-A[2])^2) >= rA &
                    sqrt((x-B[1])^2+(y-B[2])^2) >= rB)
    if (length(keep) < 2) next
    kx <- x[keep]; ky <- y[keep]
    end  <- c(kx[length(kx)], ky[length(ky)])
    prev <- c(kx[length(kx)-1], ky[length(ky)-1])
    ue <- end - prev; ue <- ue / sqrt(sum(ue^2))
    # trim the drawn line a touch before the head so it doesn't poke through
    back <- if (!is.na(s) && s > 0) headlen else 0
    nline <- max(2, length(kx) - 1)
    lines(c(kx[1:nline]) - ue[1]*back, c(ky[1:nline]) - ue[2]*back, col = col, lwd = lwd)
    heads[[length(heads)+1]] <- list(end = end, ue = ue, s = s, col = col)
  }
  
  # --- nodes + labels ---
  th <- seq(0, 2*pi, length.out = 72)
  for (i in seq_len(nrow(coords))) {
    polygon(coords[i,1] + radii[i]*cos(th), coords[i,2] + radii[i]*sin(th),
            col = fills[i], border = "grey30", lwd = 1.4)
    text(coords[i,1], coords[i,2], labels[i], font = 2, cex = label_cex[i])
  }
  
  # --- pass 2: heads on top ---
  for (h in heads) {
    end <- h$end; ue <- h$ue; s <- h$s; col <- h$col; pp <- c(-ue[2], ue[1])
    if (is.na(s) || s > 0) {                      # arrowhead (activation / unknown)
      tip  <- end
      base <- tip - ue*headlen
      polygon(c(tip[1], base[1]+pp[1]*headwid, base[1]-pp[1]*headwid),
              c(tip[2], base[2]+pp[2]*headwid, base[2]-pp[2]*headwid),
              col = col, border = col)
    } else {                                      # flat T-bar (inhibition)  ---|
      p1 <- end + pp*barhalf; p2 <- end - pp*barhalf
      segments(p1[1], p1[2], p2[1], p2[2], col = col, lwd = lwd + 1.5)
    }
  }
}

legend_signed <- function() {
  legend("bottomleft", bty = "n", cex = 0.95,
         legend = c("Requested protein", "Intermediate protein",
                    "Activation  (--->)", "Inhibition  (---|)"),
         pch = c(21, 21, NA, NA), pt.bg = c(MAIN_FILL, INTER_FILL, NA, NA),
         pt.cex = c(1.8, 1.4, NA, NA), lty = c(NA, NA, 1, 1),
         lwd = c(NA, NA, 3, 3), col = c("grey30","grey30", ACT_COL, INH_COL))
}

# =============================================================================
# (A) DETAILED
# =============================================================================
sub_edges <- unique(do.call(rbind, lapply(detailed_routes, function(vp)
  data.frame(from = vp[-length(vp)], to = vp[-1], stringsAsFactors = FALSE))))
sub_edges$sign <- sign_of(sub_edges$from, sub_edges$to)
sub <- igraph::graph_from_data_frame(sub_edges, directed = TRUE)
vn  <- igraph::V(sub)$name
is_main <- vn %in% three_nodes

set.seed(1)
co <- igraph::layout_with_fr(sub)
sp <- max(diff(range(co[,1])), diff(range(co[,2]))); sp <- ifelse(sp == 0, 1, sp)
co <- cbind((co[,1]-mean(range(co[,1])))/sp*2, (co[,2]-mean(range(co[,2])))/sp*2)

idx <- match(sub_edges$from, vn); jdx <- match(sub_edges$to, vn)
e_df <- data.frame(fi = idx, ti = jdx, sign = sub_edges$sign)
radii    <- ifelse(is_main, 0.16, 0.10)
fills    <- ifelse(is_main, MAIN_FILL, INTER_FILL)
label_cex<- ifelse(is_main, 0.95, 0.62)

draw_detailed <- function() {
  draw_signed_graph(co, vn, radii, fills, label_cex, e_df,
                    title = paste("Shortest routes among:", paste(three_nodes, collapse = ", ")),
                    curv = 0.12, headlen = 0.11, headwid = 0.06, barhalf = 0.09, lwd = 2.3, pad = 0.35)
  legend_signed()
}
png(paste0(out_detailed, ".png"), 1400, 1100, res = 150); draw_detailed(); dev.off()
pdf(paste0(out_detailed, ".pdf"), 9, 7);                   draw_detailed(); dev.off()

# =============================================================================
# (B) SIMPLIFIED NET-EFFECT TOPOLOGY
# =============================================================================
net_df <- do.call(rbind, lapply(names(chosen_routes), function(k) {
  r <- chosen_routes[[k]]; data.frame(from = r[1], to = r[length(r)],
                                      sign = net_effect(r), stringsAsFactors = FALSE) }))
net_df <- net_df[!is.na(net_df$sign), ]

ang <- pi/2 + c(0, 2*pi/3, 4*pi/3)                 # equilateral triangle
co3 <- cbind(cos(ang), sin(ang))
rownames(co3) <- three_nodes
e3  <- data.frame(fi = match(net_df$from, three_nodes),
                  ti = match(net_df$to,   three_nodes), sign = net_df$sign)

draw_simple <- function() {
  draw_signed_graph(co3, three_nodes, rep(0.34, 3), rep(MAIN_FILL, 3), rep(1.2, 3), e3,
                    title = paste("Simplified net effect among:", paste(three_nodes, collapse = ", ")),
                    curv = 0.22, headlen = 0.17, headwid = 0.10, barhalf = 0.15, lwd = 3.4, pad = 0.5)
  legend("bottomleft", bty = "n", cex = 0.95,
         legend = c("Net activation  (--->)", "Net inhibition  (---|)"),
         lty = 1, lwd = 3.4, col = c(ACT_COL, INH_COL))
}
png(paste0(out_simple, ".png"), 1100, 1000, res = 150); draw_simple(); dev.off()
pdf(paste0(out_simple, ".pdf"), 7, 6.5);                 draw_simple(); dev.off()

cat("Saved detailed:   ", out_detailed, ".png / .pdf  (",
    igraph::vcount(sub), " nodes, ", nrow(e_df), " edges)\n", sep = "")
cat("Saved simplified: ", out_simple,  ".png / .pdf  (", nrow(net_df), " net edges)\n", sep = "")