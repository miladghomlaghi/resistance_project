# =============================================================================
# plot_three_node_paths.R
#
# Given THREE proteins (nodes), using the same signed, directed interaction
# network as the topology-search pipeline, this script produces:
#
#   (A) DETAILED path plot  -> three_node_paths.png / .pdf
#       The SHORTEST route between each ordered pair, INCLUDING the intermediate
#       ("middle") proteins. Edges coloured by sign (+1 green, -1 red).
#
#   (B) SIMPLIFIED net-effect topology -> three_node_net_topology.png / .pdf
#       Each intermediate chain A -> ... -> B is collapsed to ONE edge whose
#       sign is the PRODUCT of the edge signs along the shortest route
#       (= net activation / inhibition). Also prints a 3x3 net-effect matrix
#       in the pipeline's orientation (row = from, col = to).
#
# How this maps to the pipeline's own functions:
#   - F02 get_link_topology_extract(): k_shortest_paths() sorted shortest-first
#         -> we take the shortest route (route_mode = "shortest").
#   - F19 calculate_path_value(..., full = FALSE): prod() of link signs along a
#         path -> that product is exactly the "net effect" used below.
#
# Run from the project root (beside the .RData file):
#   Rscript plot_three_node_paths.R MAPK3 AKT3 CHUK
# or edit the CONFIG block and: source("plot_three_node_paths.R")
# =============================================================================

# ---------------------------- CONFIG -----------------------------------------
three_nodes <- c("MAPK3", "AKT3", "CHUK")   # default example; override via args
data_file   <- "new_initial_info_9_6_phospho_bind.RData"
out_detailed <- "three_node_paths"           # detailed plot prefix
out_simple   <- "three_node_net_topology"     # simplified plot prefix

# How to choose the route for each pair:
#   "shortest" -> the single shortest route (mirrors F02 shortest-first). [default]
#   "all"      -> keep ALL shortest routes in the DETAILED plot.
# The SIMPLIFIED topology always uses the single shortest route (well-defined).
route_mode <- "shortest"
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) three_nodes <- args[1:3]
if (length(three_nodes) != 3) stop("Please supply exactly three protein names.")

if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("The 'igraph' package is required. Install it with install.packages('igraph').")
}
library(igraph)

# ---- Load the network -------------------------------------------------------
if (!file.exists(data_file)) stop("Cannot find '", data_file, "'. Run from project root.")
load(data_file)                      # provides: nodes, links, targets, net
if (!exists("links")) stop("'links' not found in ", data_file)

# Read links POSITIONALLY (as the pipeline does): 1=from, 2=to, 3=sign.
edges <- data.frame(
  from = as.character(links[[1]]),
  to   = as.character(links[[2]]),
  sign = suppressWarnings(as.numeric(links[[3]])),
  stringsAsFactors = FALSE
)
edges <- edges[stats::complete.cases(edges[, c("from", "to", "sign")]), ]
edges <- edges[edges$from != edges$to, ]

g <- igraph::graph_from_data_frame(edges, directed = TRUE)   # full signed graph

# Fast lookup: sign of a directed edge from -> to
edge_key  <- paste(edges$from, edges$to, sep = "\r")
sign_of   <- function(a, b) edges$sign[match(paste(a, b, sep = "\r"), edge_key)]

# ---- Validate the three requested proteins ----------------------------------
present <- three_nodes %in% igraph::V(g)$name
if (!all(present)) {
  stop("Not in the network: ", paste(three_nodes[!present], collapse = ", "))
}

# ---- Routes between every ordered pair --------------------------------------
# shortest_route(): single shortest route (vertex names) or NULL if unreachable.
shortest_route <- function(g, a, b) {
  sp <- suppressWarnings(igraph::shortest_paths(g, from = a, to = b, mode = "out")$vpath[[1]])
  if (length(sp) < 2) return(NULL)
  igraph::V(g)$name[as.integer(sp)]
}
# all_shortest_routes(): list of all shortest routes.
all_shortest_routes <- function(g, a, b) {
  res <- suppressWarnings(igraph::all_shortest_paths(g, from = a, to = b, mode = "out")$res)
  res <- Filter(function(p) length(p) > 1, res)
  lapply(res, function(p) igraph::V(g)$name[as.integer(p)])
}
# net effect of a route = product of edge signs along it (cf. F19, full=FALSE)
net_effect <- function(route) {
  s <- sign_of(route[-length(route)], route[-1])
  if (any(is.na(s))) return(NA_real_)
  prod(s)
}

pair_grid <- expand.grid(a = three_nodes, b = three_nodes, stringsAsFactors = FALSE)
pair_grid <- pair_grid[pair_grid$a != pair_grid$b, ]

chosen_routes   <- list()   # one shortest route per reachable pair (for net topology)
detailed_routes <- list()   # routes to draw in the detailed plot
net_matrix <- matrix(NA_real_, 3, 3, dimnames = list(three_nodes, three_nodes))
report <- character(0)

for (i in seq_len(nrow(pair_grid))) {
  a <- pair_grid$a[i]; b <- pair_grid$b[i]
  one <- shortest_route(g, a, b)
  if (is.null(one)) {
    report <- c(report, sprintf("  %s -> %s : no directed path", a, b))
    next
  }
  ne <- net_effect(one)
  net_matrix[a, b] <- ne
  chosen_routes[[paste(a, b)]] <- one
  mid <- if (length(one) > 2) paste(one[-c(1, length(one))], collapse = " -> ") else "(direct)"
  report <- c(report, sprintf("  %s -> %s : %d hop(s), net = %s | middle: %s",
                              a, b, length(one) - 1,
                              ifelse(is.na(ne), "NA", ifelse(ne > 0, "+ (activation)", "- (inhibition)")),
                              mid))
  if (route_mode == "all") {
    detailed_routes <- c(detailed_routes, all_shortest_routes(g, a, b))
  } else {
    detailed_routes <- c(detailed_routes, list(one))
  }
}

cat("Route summary (shortest route per ordered pair):\n",
    paste(report, collapse = "\n"), "\n\n", sep = "")
cat("Net-effect matrix (row = from, col = to; +1 activation, -1 inhibition, NA none):\n")
print(net_matrix)
cat("\n")

if (length(detailed_routes) == 0) stop("No directed routes between any pair.")

# =============================================================================
# (A) DETAILED PLOT -- routes with intermediate proteins
# =============================================================================
sub_edges <- unique(do.call(rbind, lapply(detailed_routes, function(vp)
  data.frame(from = vp[-length(vp)], to = vp[-1], stringsAsFactors = FALSE))))
sub_edges$sign <- sign_of(sub_edges$from, sub_edges$to)
sub <- igraph::graph_from_data_frame(sub_edges, directed = TRUE)

vn <- igraph::V(sub)$name
is_main <- vn %in% three_nodes
igraph::V(sub)$color       <- ifelse(is_main, "#E8833A", "#9ECAE1")
igraph::V(sub)$size        <- ifelse(is_main, 26, 15)
igraph::V(sub)$label.cex   <- ifelse(is_main, 1.0, 0.72)
igraph::V(sub)$label.font  <- ifelse(is_main, 2, 1)
igraph::V(sub)$frame.color <- "grey30"
igraph::V(sub)$label.color <- "black"
igraph::E(sub)$color <- ifelse(is.na(igraph::E(sub)$sign), "grey60",
                               ifelse(igraph::E(sub)$sign > 0, "#2CA25F", "#DE2D26"))
igraph::E(sub)$width <- 2; igraph::E(sub)$arrow.size <- 0.5

set.seed(1); lay <- igraph::layout_with_fr(sub)
draw_detailed <- function() {
  plot(sub, layout = lay, vertex.label = vn, edge.curved = 0.1,
       main = paste("Shortest routes among:", paste(three_nodes, collapse = ", ")))
  legend("bottomleft", bty = "n", cex = 0.9,
         legend = c("Requested protein", "Intermediate protein",
                    "Activation (+1)", "Inhibition (-1)"),
         pch = c(21, 21, NA, NA), pt.bg = c("#E8833A", "#9ECAE1", NA, NA),
         pt.cex = c(1.8, 1.4, NA, NA), lty = c(NA, NA, 1, 1),
         lwd = c(NA, NA, 3, 3), col = c("grey30", "grey30", "#2CA25F", "#DE2D26"))
}
png(paste0(out_detailed, ".png"), 1400, 1100, res = 150); draw_detailed(); dev.off()
pdf(paste0(out_detailed, ".pdf"), 9, 7);                   draw_detailed(); dev.off()

# =============================================================================
# (B) SIMPLIFIED NET-EFFECT TOPOLOGY -- intermediates collapsed to one edge
# =============================================================================
net_df <- do.call(rbind, lapply(names(chosen_routes), function(k) {
  r <- chosen_routes[[k]]; ne <- net_effect(r)
  data.frame(from = r[1], to = r[length(r)], sign = ne, hops = length(r) - 1,
             stringsAsFactors = FALSE)
}))
net_df <- net_df[!is.na(net_df$sign), ]

simp <- igraph::graph_from_data_frame(net_df, directed = TRUE,
                                      vertices = data.frame(name = three_nodes))
igraph::V(simp)$color       <- "#E8833A"
igraph::V(simp)$size        <- 40
igraph::V(simp)$label.cex   <- 1.1
igraph::V(simp)$label.font  <- 2
igraph::V(simp)$frame.color <- "grey30"
igraph::V(simp)$label.color <- "black"
igraph::E(simp)$color      <- ifelse(igraph::E(simp)$sign > 0, "#2CA25F", "#DE2D26")
igraph::E(simp)$lty        <- 1
igraph::E(simp)$width      <- 3
igraph::E(simp)$arrow.size <- 0.8
# arrow shape: pointed for activation, "T-bar"-like (flat) for inhibition
igraph::E(simp)$arrow.mode <- 2
igraph::E(simp)$label      <- ifelse(igraph::E(simp)$sign > 0, "+", "-")
igraph::E(simp)$label.cex  <- 1.4
igraph::E(simp)$label.color <- igraph::E(simp)$color

lay2 <- igraph::layout_in_circle(simp)
draw_simple <- function() {
  plot(simp, layout = lay2, edge.curved = 0.25,
       main = paste("Simplified net effect among:", paste(three_nodes, collapse = ", ")))
  legend("bottomleft", bty = "n", cex = 0.9,
         legend = c("Net activation (product +)", "Net inhibition (product -)"),
         lty = 1, lwd = 3, col = c("#2CA25F", "#DE2D26"))
}
png(paste0(out_simple, ".png"), 1100, 1000, res = 150); draw_simple(); dev.off()
pdf(paste0(out_simple, ".pdf"), 7, 6.5);                 draw_simple(); dev.off()

cat("Saved detailed:   ", out_detailed, ".png / .pdf  (",
    igraph::vcount(sub), " nodes, ", igraph::ecount(sub), " edges)\n", sep = "")
cat("Saved simplified: ", out_simple,  ".png / .pdf  (",
    nrow(net_df), " net edges)\n", sep = "")