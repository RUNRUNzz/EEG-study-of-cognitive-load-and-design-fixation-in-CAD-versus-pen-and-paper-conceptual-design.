# NASA-TLX: CAD vs pen-and-paper
# Paired t-tests for the six NASA-TLX dimensions with Bonferroni correction.

library(tidyverse)
library(effectsize)

# Paths -------------------------------------------------------------------
data_file <- file.path("data", "NASA_semantic.csv")
results_dir <- "results"
figures_dir <- "figures"
dir.create(results_dir, showWarnings = FALSE)
dir.create(figures_dir, showWarnings = FALSE)

# Load and clean -----------------------------------------------------------
nasa <- read_csv(data_file, show_col_types = FALSE) %>%
  mutate(
    Participant = trimws(as.character(Participant)),
    Condition = toupper(trimws(as.character(Tool))),
    Condition = recode(Condition, "P&P" = "PP"),
    across(starts_with("NASA_"), as.numeric)
  ) %>%
  filter(Condition %in% c("CAD", "PP"))

nasa_dims <- paste0("NASA_", 1:6)
nasa_labels <- c(
  NASA_1 = "Mental Demand",
  NASA_2 = "Physical Demand",
  NASA_3 = "Temporal Demand",
  NASA_4 = "Performance",
  NASA_5 = "Effort",
  NASA_6 = "Frustration"
)

required <- c("Participant", "Condition", nasa_dims)
missing <- setdiff(required, names(nasa))
if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

# One row per participant and condition ----------------------------------
nasa <- nasa %>%
  group_by(Participant, Condition) %>%
  summarise(across(all_of(nasa_dims), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# Paired tests -------------------------------------------------------------
run_nasa_test <- function(var) {
  paired <- nasa %>%
    select(Participant, Condition, Score = all_of(var)) %>%
    pivot_wider(names_from = Condition, values_from = Score) %>%
    drop_na(CAD, PP)

  test <- t.test(paired$CAD, paired$PP, paired = TRUE)
  d <- effectsize::cohens_d(paired$CAD, paired$PP, paired = TRUE, ci = 0.95)

  tibble(
    variable = var,
    Dimension = nasa_labels[[var]],
    n = nrow(paired),
    CAD_M = mean(paired$CAD),
    CAD_SD = sd(paired$CAD),
    PP_M = mean(paired$PP),
    PP_SD = sd(paired$PP),
    Mean_Difference = mean(paired$CAD - paired$PP),
    t = unname(test$statistic),
    df = unname(test$parameter),
    p = test$p.value,
    Cohens_dz = d$Cohens_d,
    dz_CI_low = d$CI_low,
    dz_CI_high = d$CI_high
  )
}

nasa_results <- map_dfr(nasa_dims, run_nasa_test) %>%
  mutate(p_adj = p.adjust(p, method = "bonferroni"))

write_csv(nasa_results, file.path(results_dir, "NASA_TLX_results.csv"))
print(nasa_results)

# Figure ------------------------------------------------------------------
plot_data <- map_dfr(nasa_dims, function(var) {
  nasa %>%
    select(Participant, Condition, Score = all_of(var)) %>%
    pivot_wider(names_from = Condition, values_from = Score) %>%
    drop_na(CAD, PP) %>%
    pivot_longer(c(CAD, PP), names_to = "Condition", values_to = "Score") %>%
    mutate(variable = var)
}) %>%
  mutate(
    Dimension = recode(variable, !!!nasa_labels),
    Dimension = factor(Dimension, levels = unname(nasa_labels)),
    Condition = factor(Condition, levels = c("CAD", "PP"),
                       labels = c("CAD", "Pen-and-paper"))
  )

plot_summary <- plot_data %>%
  group_by(variable, Dimension, Condition) %>%
  summarise(
    n = n(),
    mean = mean(Score),
    se = sd(Score) / sqrt(n),
    .groups = "drop"
  )

p_labels <- plot_summary %>%
  group_by(Dimension) %>%
  summarise(y = max(mean + se) + 2.2, .groups = "drop") %>%
  left_join(
    nasa_results %>%
      transmute(
        Dimension,
        label = if_else(p_adj < .001, "p_adj < .001",
                        paste0("p_adj = ", sprintf("%.3f", p_adj)))
      ),
    by = "Dimension"
  )

condition_colors <- c("CAD" = "#4C78A8", "Pen-and-paper" = "#F58518")
condition_shapes <- c("CAD" = 21, "Pen-and-paper" = 24)
pd <- position_dodge(width = 0.72)

nasa_plot <- ggplot(plot_summary, aes(Dimension, mean, fill = Condition)) +
  geom_col(position = pd, width = 0.62, color = "black", linewidth = 0.45) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se, group = Condition),
    position = pd, width = 0.16, linewidth = 0.65
  ) +
  geom_point(
    aes(shape = Condition, group = Condition),
    position = pd, size = 3.5, fill = "white", color = "black"
  ) +
  geom_text(
    data = p_labels,
    aes(Dimension, y, label = label),
    inherit.aes = FALSE, size = 4
  ) +
  scale_fill_manual(values = condition_colors, name = "Design medium") +
  scale_shape_manual(values = condition_shapes, name = "Design medium") +
  scale_y_continuous(
    limits = c(0, max(p_labels$y) + 1.2),
    breaks = seq(0, 20, 5),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = NULL, y = "NASA-TLX rating") +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 35, hjust = 1, size = 12.5),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4)
  )

print(nasa_plot)

ggsave(file.path(figures_dir, "NASA_TLX_CAD_PP.png"), nasa_plot,
       width = 10.5, height = 5.8, dpi = 600, bg = "white")
ggsave(file.path(figures_dir, "NASA_TLX_CAD_PP.pdf"), nasa_plot,
       width = 10.5, height = 5.8)
