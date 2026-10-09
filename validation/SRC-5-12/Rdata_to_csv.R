# Set your folder paths
folder      <- "../../Results/length_7_phospho_bind_7_6/"
csv_folder  <- "../../Results/length_7_phospho_bind_7_6/csv/"
txt_folder  <- "../../Results/length_7_phospho_bind_7_6/txt/"

# Make sure subfolders exist
if (!dir.exists(csv_folder)) dir.create(csv_folder, recursive = TRUE)
if (!dir.exists(txt_folder)) dir.create(txt_folder, recursive = TRUE)

# List all .RData files in the folder
rdata_files <- list.files(folder, pattern = "\\.RData$", full.names = TRUE)

list_csv <- character()


all_rebounded_list_numbers = data.frame(target=NA,`node3` =NA,number_of_node2=NA)



for (f in rdata_files) {
  # Load object name from RData (returns the name(s) of loaded objects)
  obj_name <- load(f)
  obj_name <- obj_name[1]        # if multiple, take first
  
  # Grab that object from the environment
  dat <- get(obj_name)
  
  # Many of your files seem to have a list with a data frame in [[1]]
  # Adjust here depending on structure:
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
  
  
  # Base name for both csv and txt
  base_name <- tools::file_path_sans_ext(basename(f))
  list_csv  <- c(list_csv, base_name)
  # Paths
  csv_path <- file.path(csv_folder, paste0(base_name, ".csv"))
  txt_path <- file.path(txt_folder,  paste0(base_name, "_col3.txt"))
  
  # Write full data frame
  write.csv(df, csv_path, row.names = FALSE)
  
  # Write third column values (one per line)
  write.csv(
    rebounded_list,
    txt_path,
    row.names = FALSE,

  )

}

# ---- pick top 40 y items (node3) overall ----
top_y_nodes <- all_rebounded_list_numbers %>%
  group_by(node3) %>%
  summarise(total_freq = sum(number_of_node2, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_freq)) %>%
  slice_head(n = 40) %>%
  pull(node3)

# ---- pick top 50 x items (target) overall ----
top_x_targets <- all_rebounded_list_numbers %>%
  group_by(target) %>%
  summarise(total_freq = sum(number_of_node2, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_freq)) %>%
  slice_head(n = 50) %>%
  pull(target)

# ---- filter to top 40 y + top 50 x ----
df_filt <- all_rebounded_list_numbers %>%
  filter(
    node3 %in% top_y_nodes,
    target %in% top_x_targets
  )

# ---- order axes by total frequency (for readability) ----
y_order <- df_filt %>%
  group_by(node3) %>%
  summarise(total = sum(number_of_node2, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(node3)

x_order <- df_filt %>%
  group_by(target) %>%
  summarise(total = sum(number_of_node2, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(target)

df_filt <- df_filt %>%
  mutate(
    node3  = factor(node3, levels = y_order),
    target = factor(target, levels = x_order)
  )

# ---- bubble heatmap ----
p <- ggplot(
  df_filt,
  aes(
    x = target,
    y = node3,
    size = number_of_node2,
    colour = number_of_node2
  )
) +
  geom_point(alpha = 0.9) +
  scale_size_continuous(
    range = c(0.4, 3.5),   # smaller circles
    trans = "sqrt"         # compress big values
  ) +
  scale_colour_viridis_c(
    option = "magma",
    direction = -1
  ) +
  labs(
    x = "Target",
    y = "Node 3",
    size = "Frequency",
    colour = "Frequency"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title  = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    panel.grid.major = element_line(colour = "grey90"),
    panel.grid.minor = element_blank()
  )

p

ggsave(
  filename = "target_node3_bubble_heatmap.pdf",
  plot     = p,
  device   = cairo_pdf,   # better text rendering
  width    = 16,          # adjust if needed
  height   = 12,
  units    = "in"
)

all_rebounded_topos = topos[frec$Topology,]
all_rebounded_topos$Topology = frec$Topology
all_rebounded_topos$frequency_perc = round(frec$Mean*100,2)

robust_topos = topos[frecb$Topology,]
robust_topos$Topology = frecb$Topology
robust_topos$frequency_perc = round(frecb$Mean*100,2)
write.csv(robust_topos, "all_rebounded_topos.csv", row.names = FALSE)
write.csv(all_rebounded_topos, "robust_topos.csv", row.names = FALSE)
