rm()
# DATA TRANSFORMATION AND NEW VARIABLES -----------------------------------
############################################################


library(igraph)
library(prodlim)
library(openxlsx)
library(foreach)
library(doParallel)


num_nodes <- 3






###################### creating the subnetwork based on the target protein

load("Results/length_7_phospho_bind_7_6/EGFR.RData")
output[[1]]


# install.packages(c("ggalluvial", "tidyverse"))
library(tidyverse)
library(ggalluvial)

# # --- Example data like yours (100 rows) ---
# set.seed(1)
# 
# df <- tibble(
#   col1 = sample(c("X", "Y", "Z"), 100, replace = TRUE, prob = c(0.45, 0.35, 0.20)),
#   # Make it so X links mostly to just 2 col2 names
#   col2 = case_when(
#     col1 == "X" ~ sample(c("B1","B2"), 100, replace = TRUE, prob = c(0.7, 0.3)),
#     TRUE        ~ sample(c("B1","B2","B3","B4","B5"), 100, replace = TRUE)
#   ),
#   # col3 depends somewhat on col2 (so B1/B2 tend to map to 3 col3 names)
#   col3 = case_when(
#     col2 %in% c("B1","B2") ~ sample(c("C1","C2","C3"), 100, replace = TRUE),
#     TRUE                  ~ sample(c("C1","C2","C3","C4","C5"), 100, replace = TRUE)
#   )
# )


df = output[[1]]
colnames(df)=c("Node 1","Node 2","Node 3")
# --- Aggregate to counts for plotting (recommended) ---
plot_df <- df %>%
  count(`Node 1`, `Node 2`, `Node 3`, name = "n")

p <- ggplot(plot_df[1:100,],
            aes(axis1 = `Node 1`, axis2 = `Node 2`, axis3 = `Node 3`, y = n)) +
  geom_alluvium(aes(fill = `Node 1`), alpha = 0.75, width = 0.15) +
  geom_stratum(width = 0.22, color = "grey30", fill = "grey95") +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 3) +
  scale_x_discrete(limits = c("Node 1", "Node 2", "Node 3"),
                   expand = c(0.05, 0.05)) +
  labs(x = NULL, y = "Count") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.ticks.x = element_blank()
  )

p

plot_df2 <- plot_df %>%
  group_by(`Node 1`) %>% mutate(t1 = sum(n)) %>% ungroup() %>%
  group_by(`Node 2`) %>% mutate(t2 = sum(n)) %>% ungroup() %>%
  group_by(`Node 3`) %>% mutate(t3 = sum(n)) %>% ungroup() %>%
  mutate(
    `Node 1` = fct_reorder(`Node 1`, t1, .desc = TRUE),
    `Node 2` = fct_reorder(`Node 2`, t2, .desc = TRUE),
    `Node 3` = fct_reorder(`Node 3`, t3, .desc = TRUE)
  )

ggplot(plot_df2[1:100,],
       aes(axis1 = `Node 1`, axis2 = `Node 2`, axis3 = `Node 3`, y = n)) +
  geom_alluvium(aes(fill = `Node 1`), alpha = 0.75, width = 0.15) +
  geom_stratum(width = 0.22, color = "grey30", fill = "grey95") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Node 1", "Node 2", "Node 3"), expand = c(0.05, 0.05)) +
  labs(x = NULL, y = "Count") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", axis.ticks.x = element_blank())
