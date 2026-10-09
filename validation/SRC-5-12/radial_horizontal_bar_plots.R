closeAllConnections()
rm(list = ls())

library(foreach)
library(pracma)
library(igraph)
library(doParallel)
library(openxlsx)




load("new_initial_info_9_6_phospho_bind.RData")


# Set your folder paths
folder      <- "../../Results/length_7_phospho_bind_7_6/"


# List all .RData files in the folder
rdata_files0 <- list.files(folder, pattern = "\\.RData$", full.names = TRUE)

patterns <- c("CDK4", "FGFR2", "MTOR","SRC","MAP2K1")

rdata_files = rdata_files0[grep(paste(patterns, collapse = "|"), rdata_files0)]

list_csv <- character()


all_rebounded_list_numbers = data.frame(target=NA,`node3` =NA,number_of_node2=NA)



for (f in rdata_files) {
  obj_name <- load(f)
  obj_name <- obj_name[1]        # if multiple, take first
  
  dat <- get(obj_name)

  if (is.data.frame(dat)) {
    df <- dat
  } else if (is.list(dat) && length(dat) >= 1 && is.data.frame(dat[[1]])) {
    df <- dat[[1]]
  } else {
    message("Skipping file ", basename(f),
            " because object '", obj_name,
            "' is not a data frame or list-with-data-frame.")
    next
  }
  
  # Check it has at least 3 columns
  if (ncol(df) < 3) {
    message("Skipping file ", basename(f),
            " because data frame has < 3 columns.")
    next
  }
  
  # Extract unique values from 3rd column
  rebounded_list <- (df%>%group_by(`3`)%>%mutate(target = dat[[1]][1,1], num_X2 = n())%>%select(target,`3`,num_X2)%>%distinct()%>%arrange(desc(num_X2)))
  colnames(rebounded_list) = c("target","node3","number_of_node2")
  
  all_rebounded_list_numbers = rbind(all_rebounded_list_numbers,rebounded_list )
  
  
  
}

################################
# manually selected targets
######### Test
library(tidyverse)
library(igraph)
library(ggraph)
library(scales)



top_x_targets <- c("CDK4","FGFR2","MTOR","MAP2K1","SRC")
for (i in 1:5){
  
  
  
man_selected_rebounded_list_numbers = all_rebounded_list_numbers %>%
  filter(
    target %in% top_x_targets[i],
    
    node3 %in% nodes
  )

# ---- pick top 40 y items (node3) overall ----
top_y_nodes <- man_selected_rebounded_list_numbers %>%
  group_by(node3) %>%
  summarise(total_freq = sum(number_of_node2, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_freq)) %>%
  slice_head(n = 100) %>%
  pull(node3)

# ---- pick top selected items (target) overall ----

# ---- filter to top 40 y + top 50 x ----
df_filt <- man_selected_rebounded_list_numbers %>%
  filter(
    node3 %in% top_y_nodes,
    target %in% top_x_targets[i]
  )


# ---- Parameters you can tweak ----
R_outer   <- 1.0   # radius of circumference for outer nodes
label_pad <- 0.10  # how far labels sit beyond the circumference

# 1) Place ALL targets on the circumference (equal angles)
n <- nrow(df_filt)
angles <- seq(0, 2*pi, length.out = n + 1)[-1]

targets_xy <- tibble(
  name   = df_filt$node3 ,
  x      = cos(angles) * R_outer,
  y      = sin(angles) * R_outer,
  # angle for text rotation in degrees (radial alignment)
  angle_deg_raw = atan2(y, x) * 180 / pi
) %>%
  # Flip labels on the left half so they are not upside down
  mutate(
    angle_deg = if_else(angle_deg_raw > 90 | angle_deg_raw < -90, angle_deg_raw + 180, angle_deg_raw),
    hjust     = if_else(angle_deg_raw > 90 | angle_deg_raw < -90, 1, 0),
    # push labels slightly outward along the radius
    x_lab     = x + cos(angles) * label_pad,
    y_lab     = y + sin(angles) * label_pad
  )

layout_df <- bind_rows(
  tibble(name = unique(df_filt$target), x = 0, y = 0),  # center
  targets_xy %>% transmute(name, x, y)
)

# 2) Build graph
g <- graph_from_data_frame(df_filt %>% select(target, node3, number_of_node2), directed = FALSE)

# 3) Plot: targets on circumference; labels aligned with radius line
p <- ggraph(g, layout = "manual", x = layout_df$x, y = layout_df$y) +
  geom_edge_link(
    aes(edge_width = 1, edge_alpha = 1),
    colour = "grey40",
    lineend = "round"
  ) +
  geom_node_point(
    aes(size = if_else(name == unique(df_filt$target), 2, 2)),
    colour = "grey20"
  ) +
  # Use manual label positions (x_lab/y_lab) and rotate along radius
  geom_text(
    data = targets_xy,
    aes(x = x_lab, y = y_lab, label = name, angle = angle_deg, hjust = hjust),
    size = 3
  ) +
  scale_edge_width(range = c(0.2, 1.2), guide = "none") +
  scale_edge_alpha(range = c(0.2, 0.9), guide = "none") +
  scale_size_identity() +
  coord_equal(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25), expand = FALSE) +
  theme_void()

p <- p +
  coord_equal(
    xlim = c(-1.35, 1.35),
    ylim = c(-1.35, 1.35),
    expand = FALSE,
    clip = "off"
  ) +
  theme(
    plot.margin = margin(15, 15, 15, 15)  # top, right, bottom, left (pt)
  )
print(p)

# Optional export
ggsave(paste("radial_circumference_aligned_labels",top_x_targets[i],".png"), p, width = 6, height = 6)
}

library(tidyverse)
library(scales)

# ---- example data (replace with yours) ----
df <- tribble(
  ~inhibitor, ~highly_robust, ~robust, ~not_explained,
  "SRCi",        135,            53,          0,
  "FGFR2i",      264,           121,          1,
  "MEK1/2i",      85,            58,          0,
  "mTORC1i",     100,            57,          1,
  "CDK4/6i",      58,            24,         0
)

# ---- reshape + compute totals and explained % ----
long <- df %>%
  pivot_longer(-inhibitor, names_to = "class", values_to = "n") %>%
  mutate(
    class = recode(class,
                   highly_robust  = "Highly robust",
                   robust         = "Robust",
                   not_explained  = "Not explained")
  )

summary_df <- df %>%
  mutate(
    total = highly_robust + robust + not_explained,
    explained = highly_robust + robust,
    explained_pct = explained / total
  )

# Put inhibitors in desired order (top to bottom like the example)
long <- long %>%
  left_join(summary_df %>% select(inhibitor, total), by = "inhibitor") %>%
  mutate(inhibitor = fct_reorder(inhibitor, total, .desc = TRUE))

summary_df <- summary_df %>%
  mutate(inhibitor = factor(inhibitor, levels = levels(long$inhibitor)))

# ---- plot ----
cols <- c(
  "Highly robust"  = "#0b6b34",
  "Robust"         = "#7ff2a2",
  "Not explained"  = "#d9dde3"
)

p <- ggplot(long, aes(x = inhibitor, y = n, fill = class)) +
  geom_col(width = 0.62, color = NA) +
  
  # numbers inside segments (only if segment > 0)
  geom_text(
    data = long %>% filter(n > 0),
    aes(label = n),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    color = "white"
  ) +
  
  # percent + (N) on the right
  geom_text(
    data = summary_df,
    aes(
      x = inhibitor,
      y = total * 1.06,  # pushes label a bit to the right
      label = paste0(percent(explained_pct, accuracy = 0.1), " (", total, ")"),
      color = explained_pct
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.6,
    fontface = "bold"
  ) +
  scale_color_gradient(low = "#f59e0b", high = "#16a34a", guide = "none") +
  
  coord_flip(clip = "off") +
  scale_fill_manual(values = cols) +
  
  # add right-side space for the % labels
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  
  labs(
    x = NULL, y = NULL, fill = NULL,
    title = "Option 1: Horizontal Percentage Bars",
    subtitle = "Clean, easy to compare proportions across inhibitors"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    plot.margin = margin(10, 40, 10, 10)
  )

library(svglite)

ggsave(
  filename = "horizontal_percentage_bars.svg",
  plot     = p,
  width    = 10,
  height   = 6,
  units    = "in",
  device   = svglite,
  bg       = "white"
)

dev.off()