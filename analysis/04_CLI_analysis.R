# Baseline-corrected Cognitive Load Index (CLI): CAD vs pen-and-paper

library(tidyverse)
library(effectsize)
library(patchwork)

# Paths -------------------------------------------------------------------
data_file <- file.path("data", "PSD_CLI.csv")
results_dir <- "results"
figures_dir <- "figures"
dir.create(results_dir, showWarnings = FALSE)
dir.create(figures_dir, showWarnings = FALSE)

# Load and clean -----------------------------------------------------------
eeg <- read_csv(data_file, show_col_types = FALSE) %>%
  transmute(
    Participant = trimws(as.character(Participant)),
    Condition = toupper(trimws(as.character(Condition))),
    Phase = tolower(trimws(as.character(Phase))),
    CLI = as.numeric(CLI)
  ) %>%
  filter(Participant != "", !is.na(Participant),
         Condition %in% c("CAD", "PP"),
         Phase %in% c("baseline", "free_design"))

# Check one value per participant-condition-phase -------------------------
duplicates <- eeg %>%
  count(Participant, Condition, Phase) %>%
  filter(n > 1)
if (nrow(duplicates) > 0) stop("Duplicate participant-condition-phase rows found.")

# Baseline correction ------------------------------------------------------
cli_corrected <- eeg %>%
  pivot_wider(names_from = Phase, values_from = CLI) %>%
  mutate(Delta_CLI = free_design - baseline) %>%
  select(Participant, Condition, baseline, free_design, Delta_CLI)

cli_paired <- cli_corrected %>%
  select(Participant, Condition, Delta_CLI) %>%
  pivot_wider(names_from = Condition, values_from = Delta_CLI) %>%
  drop_na(CAD, PP) %>%
  mutate(Difference_CAD_minus_PP = CAD - PP)

# Statistics ---------------------------------------------------------------
cli_descriptives <- bind_rows(
  cli_paired %>% summarise(Condition = "CAD", n = n(), M = mean(CAD), SD = sd(CAD)),
  cli_paired %>% summarise(Condition = "PP", n = n(), M = mean(PP), SD = sd(PP))
)

shapiro_result <- shapiro.test(cli_paired$Difference_CAD_minus_PP)
cli_ttest <- t.test(cli_paired$CAD, cli_paired$PP, paired = TRUE, conf.level = 0.95)
cli_dz <- effectsize::cohens_d(cli_paired$CAD, cli_paired$PP, paired = TRUE, ci = 0.95)
cli_wilcoxon <- wilcox.test(cli_paired$CAD, cli_paired$PP,
                            paired = TRUE, exact = FALSE)

cli_results <- tibble(
  Paired_n = nrow(cli_paired),
  CAD_M = mean(cli_paired$CAD),
  CAD_SD = sd(cli_paired$CAD),
  PP_M = mean(cli_paired$PP),
  PP_SD = sd(cli_paired$PP),
  Mean_Difference_CAD_minus_PP = mean(cli_paired$Difference_CAD_minus_PP),
  Difference_CI_low = cli_ttest$conf.int[1],
  Difference_CI_high = cli_ttest$conf.int[2],
  t = unname(cli_ttest$statistic),
  df = unname(cli_ttest$parameter),
  p = cli_ttest$p.value,
  Cohens_dz = cli_dz$Cohens_d,
  dz_CI_low = cli_dz$CI_low,
  dz_CI_high = cli_dz$CI_high,
  Shapiro_W = unname(shapiro_result$statistic),
  Shapiro_p = shapiro_result$p.value,
  Wilcoxon_V = unname(cli_wilcoxon$statistic),
  Wilcoxon_p = cli_wilcoxon$p.value
)

write_csv(cli_corrected, file.path(results_dir, "CLI_baseline_corrected.csv"))
write_csv(cli_paired, file.path(results_dir, "CLI_complete_pairs.csv"))
write_csv(cli_descriptives, file.path(results_dir, "CLI_descriptives.csv"))
write_csv(cli_results, file.path(results_dir, "CLI_results.csv"))

print(cli_results)

# Figure ------------------------------------------------------------------
plot_long <- cli_paired %>%
  select(Participant, CAD, PP) %>%
  pivot_longer(c(CAD, PP), names_to = "Condition", values_to = "Delta_CLI") %>%
  mutate(Condition = factor(Condition, levels = c("PP", "CAD"),
                            labels = c("Pen-and-paper", "CAD")))

condition_summary <- plot_long %>%
  group_by(Condition) %>%
  summarise(
    n = n(),
    Mean = mean(Delta_CLI),
    SD = sd(Delta_CLI),
    SE = SD / sqrt(n),
    CI_low = Mean - qt(.975, n - 1) * SE,
    CI_high = Mean + qt(.975, n - 1) * SE,
    .groups = "drop"
  )

condition_colors <- c("Pen-and-paper" = "#F58518", "CAD" = "#4C78A8")
condition_shapes <- c("Pen-and-paper" = 17, "CAD" = 16)

panel_a <- ggplot(plot_long, aes(Condition, Delta_CLI, group = Participant)) +
  geom_line(color = "grey75", linewidth = 0.5, alpha = 0.65) +
  geom_point(aes(color = Condition, shape = Condition), size = 2.8, alpha = 0.85) +
  geom_errorbar(
    data = condition_summary,
    aes(Condition, ymin = CI_low, ymax = CI_high),
    inherit.aes = FALSE, width = 0.10, linewidth = 0.85
  ) +
  geom_point(
    data = filter(condition_summary, Condition == "Pen-and-paper"),
    aes(Condition, Mean), inherit.aes = FALSE,
    shape = 24, size = 4.7, fill = "white", color = "black", stroke = 1.1
  ) +
  geom_point(
    data = filter(condition_summary, Condition == "CAD"),
    aes(Condition, Mean), inherit.aes = FALSE,
    shape = 21, size = 4.7, fill = "white", color = "black", stroke = 1.1
  ) +
  scale_color_manual(values = condition_colors) +
  scale_shape_manual(values = condition_shapes) +
  labs(x = NULL, y = expression(Delta * "CLI (Free design - Baseline)")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

diff_summary <- tibble(
  Mean = mean(cli_paired$Difference_CAD_minus_PP),
  CI_low = cli_ttest$conf.int[1],
  CI_high = cli_ttest$conf.int[2]
)

p_label <- if (cli_ttest$p.value < .001) {
  "p < .001"
} else {
  paste0("p = ", sub("^0", "", sprintf("%.3f", cli_ttest$p.value)))
}

dz_value <- mean(cli_paired$Difference_CAD_minus_PP) /
  sd(cli_paired$Difference_CAD_minus_PP)

result_label <- paste0(
  "Mean difference = ", sprintf("%.2f", diff_summary$Mean),
  "\n95% CI [", sprintf("%.2f", diff_summary$CI_low), ", ",
  sprintf("%.2f", diff_summary$CI_high), "]",
  "\n", p_label, ", dz = ", sprintf("%.2f", dz_value)
)

panel_b <- ggplot(cli_paired, aes(x = 1, y = Difference_CAD_minus_PP)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_jitter(width = 0.08, size = 2.6, alpha = 0.65, color = "grey40") +
  geom_errorbar(
    data = diff_summary,
    aes(x = 1, ymin = CI_low, ymax = CI_high),
    inherit.aes = FALSE, width = 0.10, linewidth = 1
  ) +
  geom_point(
    data = diff_summary,
    aes(x = 1, y = Mean),
    inherit.aes = FALSE, shape = 23, size = 5,
    fill = "white", color = "black", stroke = 1.1
  ) +
  annotate(
    "text", x = 1.28,
    y = max(cli_paired$Difference_CAD_minus_PP, na.rm = TRUE),
    label = result_label, hjust = 0, vjust = 1, size = 4
  ) +
  scale_x_continuous(breaks = 1, labels = "CAD - PP", limits = c(0.75, 1.75)) +
  labs(x = NULL, y = expression(Delta * "CLI difference")) +
  theme_classic(base_size = 14)

cli_plot <- panel_a + panel_b +
  plot_layout(widths = c(1.35, 1)) +
  plot_annotation(tag_levels = "A")

print(cli_plot)

ggsave(file.path(figures_dir, "CLI_CAD_PP_publication_plot.png"), cli_plot,
       width = 10, height = 5.2, dpi = 600, bg = "white")
ggsave(file.path(figures_dir, "CLI_CAD_PP_publication_plot.pdf"), cli_plot,
       width = 10, height = 5.2)
