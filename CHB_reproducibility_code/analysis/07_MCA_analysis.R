# Multiple Correspondence Analysis (MCA) of categorized EEG TRP states

library(tidyverse)
library(FactoMineR)
library(factoextra)

# Paths -------------------------------------------------------------------
data_file <- file.path("data", "MCA_ERS_ERD.csv")
results_dir <- "results"
figures_dir <- "figures"
dir.create(results_dir, showWarnings = FALSE)
dir.create(figures_dir, showWarnings = FALSE)

# Load the categorized TRP data -------------------------------------------
mca_df <- read_csv(data_file, show_col_types = FALSE) %>%
  mutate(
    Condition = toupper(trimws(as.character(Condition))),
    Condition = recode(Condition, "P&P" = "PP"),
    Condition = factor(Condition, levels = c("CAD", "PP"))
  ) %>%
  filter(!is.na(Condition))

state_cols <- setdiff(names(mca_df), "Condition")
mca_df <- mca_df %>%
  mutate(across(all_of(state_cols), factor))

# Remove any feature that has fewer than two observed states.
keep_cols <- state_cols[sapply(mca_df[state_cols], function(x) {
  nlevels(droplevels(x)) >= 2
})]
mca_df <- mca_df %>% select(Condition, all_of(keep_cols))

# MCA: Condition is supplementary ----------------------------------------
set.seed(42)
res_mca <- MCA(mca_df, quali.sup = 1, graph = FALSE)

# Save core results -------------------------------------------------------
eigenvalues <- as.data.frame(res_mca$eig) %>%
  rownames_to_column("Dimension")
condition_coordinates <- as.data.frame(res_mca$quali.sup$coord) %>%
  rownames_to_column("Condition")
condition_vtests <- as.data.frame(res_mca$quali.sup$v.test) %>%
  rownames_to_column("Condition")
category_contributions <- as.data.frame(res_mca$var$contrib) %>%
  rownames_to_column("Category")

write_csv(eigenvalues, file.path(results_dir, "MCA_eigenvalues.csv"))
write_csv(condition_coordinates, file.path(results_dir, "MCA_condition_coordinates.csv"))
write_csv(condition_vtests, file.path(results_dir, "MCA_condition_vtests.csv"))
write_csv(category_contributions, file.path(results_dir, "MCA_category_contributions.csv"))

cat("Dimension 1:", round(res_mca$eig[1, 2], 2), "%\n")
cat("Dimension 2:", round(res_mca$eig[2, 2], 2), "%\n")
print(round(res_mca$quali.sup$coord[, 1:2], 3))
print(round(res_mca$quali.sup$v.test[, 1:2], 3))

# Feature-level chi-square follow-ups ------------------------------------
cramers_v <- function(tab, chi_square) {
  n <- sum(tab)
  k <- min(nrow(tab) - 1, ncol(tab) - 1)
  if (k <= 0) return(NA_real_)
  sqrt(as.numeric(chi_square) / (n * k))
}

chi_square_results <- map_dfr(keep_cols, function(var) {
  tab <- table(mca_df[[var]], mca_df$Condition)
  test <- suppressWarnings(chisq.test(tab, correct = FALSE))

  tibble(
    variable = var,
    chi_square = unname(test$statistic),
    df = unname(test$parameter),
    p = test$p.value,
    cramers_v = cramers_v(tab, test$statistic)
  )
}) %>%
  mutate(p_adj = p.adjust(p, method = "bonferroni")) %>%
  arrange(p_adj)

write_csv(chi_square_results,
          file.path(results_dir, "MCA_feature_chi_square_results.csv"))
print(chi_square_results)

# Figure 4: MCA individual factor map ------------------------------------
dim1 <- res_mca$eig[1, 2]
dim2 <- res_mca$eig[2, 2]

coordinates <- as.data.frame(res_mca$ind$coord[, 1:2])
names(coordinates) <- c("Dim1", "Dim2")
coordinates$Condition <- factor(
  mca_df$Condition,
  levels = c("CAD", "PP"),
  labels = c("CAD", "Pen-and-paper")
)

mca_base <- fviz_mca_ind(
  res_mca,
  habillage = mca_df$Condition,
  palette = c("CAD" = "#4C78A8", "PP" = "#F58518"),
  addEllipses = TRUE,
  ellipse.type = "confidence",
  ellipse.level = 0.95,
  geom = "blank",
  title = NULL,
  legend.title = NULL
)

mca_plot <- mca_base +
  geom_point(
    data = coordinates,
    aes(Dim1, Dim2, color = Condition, shape = Condition),
    inherit.aes = FALSE,
    size = 3,
    alpha = 0.78
  ) +
  scale_color_manual(
    values = c("CAD" = "#4C78A8", "Pen-and-paper" = "#F58518"),
    name = "Design medium"
  ) +
  scale_shape_manual(
    values = c("CAD" = 16, "Pen-and-paper" = 17),
    name = "Design medium"
  ) +
  labs(
    x = paste0("Dimension 1 (", sprintf("%.1f", dim1), "%)"),
    y = paste0("Dimension 2 (", sprintf("%.1f", dim2), "%)")
  ) +
  theme_classic(base_size = 15) +
  theme(legend.position = "top")

print(mca_plot)

ggsave(file.path(figures_dir, "MCA_Map.png"), mca_plot,
       width = 8, height = 6, dpi = 600, bg = "white")
ggsave(file.path(figures_dir, "MCA_Map.pdf"), mca_plot,
       width = 8, height = 6)

# Figure 5: top 15 contributors to MCA Dimension 1 -----------------------
# Convert FactoMineR category names into the label
format_category <- function(x) {
  state <- case_when(
    str_detect(x, "ERS$") ~ "ERS",
    str_detect(x, "ERD$") ~ "ERD",
    str_detect(x, "NoChange$") ~ "NoChange",
    TRUE ~ NA_character_
  )

  feature <- x %>%
    str_remove("[._](ERS|ERD|NoChange)$") %>%
    str_remove("^TRP_") %>%
    str_remove("_diff_FD$") %>%
    str_replace_all("_", " ") %>%
    str_replace_all("Pz", "PZ") %>%
    str_replace_all("Fz", "Fz")

  paste0(feature, " (", state, ")")
}

top_contrib <- category_contributions %>%
  arrange(desc(`Dim 1`)) %>%
  slice_head(n = 15) %>%
  mutate(
    State = case_when(
      str_detect(Category, "ERS$") ~ "ERS",
      str_detect(Category, "ERD$") ~ "ERD",
      str_detect(Category, "NoChange$") ~ "NoChange",
      TRUE ~ "Other"
    ),
    Label = map_chr(Category, format_category),
    Label = factor(Label, levels = rev(Label))
  )

state_colors <- c(
  "ERD" = "#E67E52",
  "ERS" = "#5B9BD5",
  "NoChange" = "#F2C94C"
)
state_labels <- c(
  "ERD" = "ERD (Desynchronization)",
  "ERS" = "ERS (Synchronization)",
  "NoChange" = "No Change"
)

contribution_plot <- ggplot(top_contrib, aes(`Dim 1`, Label, fill = State)) +
  geom_col(width = 0.78) +
  geom_text(
    aes(label = sprintf("%.1f%%", `Dim 1`)),
    hjust = -0.12,
    size = 3.6
  ) +
  scale_fill_manual(
    values = state_colors,
    breaks = c("ERD", "ERS", "NoChange"),
    labels = state_labels,
    name = "Neural State"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(x = "Contribution (%)", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linetype = "dashed", color = "grey75"),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(face = "bold")
  )

print(contribution_plot)

ggsave(file.path(figures_dir, "MCA_Dimension1_contributions.png"), contribution_plot,
       width = 8.2, height = 5.8, dpi = 600, bg = "white")
ggsave(file.path(figures_dir, "MCA_Dimension1_contributions.pdf"), contribution_plot,
       width = 8.2, height = 5.8)
