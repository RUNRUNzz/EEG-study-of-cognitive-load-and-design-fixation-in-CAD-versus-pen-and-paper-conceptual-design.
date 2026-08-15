# Multiple Correspondence Analysis of categorized EEG TRP states

library(tidyverse)
library(FactoMineR)
library(factoextra)

# Paths -------------------------------------------------------------------
data_file <- file.path("results", "EEG_TRP_free_design_wide.csv")
results_dir <- "results"
figures_dir <- "figures"
dir.create(results_dir, showWarnings = FALSE)
dir.create(figures_dir, showWarnings = FALSE)

# Load data ---------------------------------------------------------------
df <- read_csv(data_file, show_col_types = FALSE) %>%
  mutate(Condition = factor(Condition, levels = c("CAD", "PP")))

trp_cols <- grep("^TRP_.*_diff_FD$", names(df), value = TRUE)
if (length(trp_cols) == 0) stop("No free-design TRP columns were found.")

# Categorize TRP as ERD / NoChange / ERS ---------------------------------
# This line preserves the threshold implementation used in the uploaded MCA
# script so the cleaning step does not silently change the current result.
# See REVIEW_BEFORE_RELEASE.md before the repository is finalized.
bin_trp <- function(x, k = 0.30) {
  threshold <- k * (mad(x, na.rm = TRUE) / 0.6745)
  cut(
    x,
    breaks = c(-Inf, -threshold, threshold, Inf),
    labels = c("ERD", "NoChange", "ERS"),
    include.lowest = TRUE
  )
}

df_binned <- df %>%
  mutate(across(all_of(trp_cols), bin_trp))

# Remove features with fewer than two observed states ---------------------
keep_trp <- trp_cols[sapply(df_binned[trp_cols], function(x) {
  nlevels(droplevels(factor(x))) >= 2
})]

mca_df <- df_binned %>%
  select(Condition, all_of(keep_trp)) %>%
  drop_na(Condition)

# MCA: Condition is supplementary ----------------------------------------
res_mca <- MCA(mca_df, quali.sup = 1, graph = FALSE)

# Main MCA outputs --------------------------------------------------------
eigenvalues <- as.data.frame(res_mca$eig) %>%
  rownames_to_column("Dimension")

condition_coordinates <- as.data.frame(res_mca$quali.sup$coord) %>%
  rownames_to_column("Condition")
condition_vtests <- as.data.frame(res_mca$quali.sup$v.test) %>%
  rownames_to_column("Condition")
condition_cos2 <- as.data.frame(res_mca$quali.sup$cos2) %>%
  rownames_to_column("Condition")

category_contributions <- as.data.frame(res_mca$var$contrib) %>%
  rownames_to_column("Category")
category_cos2 <- as.data.frame(res_mca$var$cos2) %>%
  rownames_to_column("Category")

write_csv(mca_df, file.path(results_dir, "MCA_binned_TRP.csv"))
write_csv(eigenvalues, file.path(results_dir, "MCA_eigenvalues.csv"))
write_csv(condition_coordinates, file.path(results_dir, "MCA_condition_coordinates.csv"))
write_csv(condition_vtests, file.path(results_dir, "MCA_condition_vtests.csv"))
write_csv(condition_cos2, file.path(results_dir, "MCA_condition_cos2.csv"))
write_csv(category_contributions, file.path(results_dir, "MCA_category_contributions.csv"))
write_csv(category_cos2, file.path(results_dir, "MCA_category_cos2.csv"))

cat("Dimension 1:", round(res_mca$eig[1, 2], 2), "%\n")
cat("Dimension 2:", round(res_mca$eig[2, 2], 2), "%\n")
print(round(res_mca$quali.sup$v.test[, 1:2], 3))

# Exploratory feature-level chi-square tests ------------------------------
# Bonferroni correction is applied across the tested TRP-state features.
cramers_v <- function(tab, chi_square) {
  n <- sum(tab)
  k <- min(nrow(tab) - 1, ncol(tab) - 1)
  if (k <= 0) return(NA_real_)
  sqrt(as.numeric(chi_square) / (n * k))
}

chi_square_results <- map_dfr(keep_trp, function(var) {
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

# Figure 4: individual MCA map --------------------------------------------
dim1 <- res_mca$eig[1, 2]
dim2 <- res_mca$eig[2, 2]

mca_coordinates <- as.data.frame(res_mca$ind$coord[, 1:2])
names(mca_coordinates) <- c("Dim1", "Dim2")
mca_coordinates$Condition <- factor(
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
    data = mca_coordinates,
    aes(Dim1, Dim2, color = Condition, shape = Condition),
    inherit.aes = FALSE, size = 3, alpha = 0.78
  ) +
  scale_color_manual(values = c("CAD" = "#4C78A8", "Pen-and-paper" = "#F58518"),
                     name = "Design medium") +
  scale_shape_manual(values = c("CAD" = 16, "Pen-and-paper" = 17),
                     name = "Design medium") +
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

# Figure 5: strongest Dimension 1 category contributions -----------------
top_contrib <- category_contributions %>%
  arrange(desc(`Dim 1`)) %>%
  slice_head(n = 10) %>%
  mutate(Category = fct_reorder(Category, `Dim 1`))

contribution_plot <- ggplot(top_contrib, aes(`Dim 1`, Category)) +
  geom_col() +
  labs(x = "Contribution to Dimension 1 (%)", y = NULL) +
  theme_classic(base_size = 14)

print(contribution_plot)

ggsave(file.path(figures_dir, "MCA_Dimension1_contributions.png"), contribution_plot,
       width = 8, height = 5.5, dpi = 600, bg = "white")
ggsave(file.path(figures_dir, "MCA_Dimension1_contributions.pdf"), contribution_plot,
       width = 8, height = 5.5)
