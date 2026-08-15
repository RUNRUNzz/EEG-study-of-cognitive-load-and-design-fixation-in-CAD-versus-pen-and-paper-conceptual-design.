# Semantic distance: Tool x furniture task mixed-effects model

library(tidyverse)
library(lmerTest)
library(emmeans)

# Paths -------------------------------------------------------------------
data_file <- file.path("data", "semantic_distance.csv")
results_dir <- "results"
figures_dir <- "figures"
dir.create(results_dir, showWarnings = FALSE)
dir.create(figures_dir, showWarnings = FALSE)

# Load and clean -----------------------------------------------------------
df <- read_csv(data_file, show_col_types = FALSE) %>%
  transmute(
    Participant = factor(trimws(as.character(Participant))),
    Tool = toupper(trimws(as.character(Tool))),
    Task = str_to_title(trimws(as.character(Task))),
    Semantic = as.numeric(Semantic)
  ) %>%
  mutate(
    Tool = factor(Tool, levels = c("CAD", "PP")),
    Task = factor(Task, levels = c("Bed", "Cabinet"))
  ) %>%
  drop_na(Participant, Tool, Task, Semantic)

# The public derived dataset should already contain the final analytic sample.
# No original participant IDs are hard-coded in this script.
if (any(count(df, Participant, Tool)$n > 1)) {
  stop("More than one participant-session row was found for a tool condition.")
}

cat("Participants:", n_distinct(df$Participant), "\n")
cat("Observations:", nrow(df), "\n")

# Sum-to-zero contrasts for Type III tests --------------------------------
options(contrasts = c("contr.sum", "contr.poly"))

# Descriptives -------------------------------------------------------------
semantic_descriptives <- df %>%
  group_by(Tool, Task) %>%
  summarise(n = n(), M = mean(Semantic), SD = sd(Semantic), .groups = "drop")

# Mixed-effects model ------------------------------------------------------
model_semantic <- lmer(
  Semantic ~ Tool * Task + (1 | Participant),
  data = df,
  REML = TRUE
)

semantic_anova <- anova(model_semantic, type = 3, ddf = "Satterthwaite")

# Estimated marginal means and planned contrasts --------------------------
tool_emm <- emmeans(model_semantic, ~ Tool, weights = "equal",
                    lmer.df = "satterthwaite")
tool_contrast <- contrast(tool_emm, list("CAD - PP" = c(1, -1)), adjust = "none")

task_emm <- emmeans(model_semantic, ~ Task, weights = "equal",
                    lmer.df = "satterthwaite")
task_contrast <- contrast(task_emm, list("Cabinet - Bed" = c(-1, 1)), adjust = "none")

cell_emm <- emmeans(model_semantic, ~ Tool * Task, lmer.df = "satterthwaite")

tool_effect <- eff_size(
  tool_emm,
  sigma = sigma(model_semantic),
  edf = df.residual(model_semantic),
  method = "pairwise"
)

# Save results -------------------------------------------------------------
write_csv(semantic_descriptives,
          file.path(results_dir, "Semantic_descriptives.csv"))
write_csv(
  as.data.frame(semantic_anova) %>% rownames_to_column("Effect"),
  file.path(results_dir, "Semantic_Type3_tests.csv")
)
write_csv(as.data.frame(summary(tool_emm, infer = c(TRUE, TRUE))),
          file.path(results_dir, "Semantic_Tool_EMMs.csv"))
write_csv(as.data.frame(summary(tool_contrast, infer = c(TRUE, TRUE))),
          file.path(results_dir, "Semantic_CAD_minus_PP.csv"))
write_csv(as.data.frame(summary(task_emm, infer = c(TRUE, TRUE))),
          file.path(results_dir, "Semantic_Task_EMMs.csv"))
write_csv(as.data.frame(summary(task_contrast, infer = c(TRUE, TRUE))),
          file.path(results_dir, "Semantic_Cabinet_minus_Bed.csv"))
write_csv(as.data.frame(summary(tool_effect, infer = c(TRUE, TRUE))),
          file.path(results_dir, "Semantic_Tool_effect_size.csv"))

print(semantic_anova)
print(summary(tool_contrast, infer = c(TRUE, TRUE)))
print(summary(task_contrast, infer = c(TRUE, TRUE)))

# Basic model diagnostics --------------------------------------------------
png(file.path(figures_dir, "Semantic_model_diagnostics.png"),
    width = 1600, height = 800, res = 180)
par(mfrow = c(1, 2))
plot(fitted(model_semantic), residuals(model_semantic),
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, lty = 2)
qqnorm(residuals(model_semantic))
qqline(residuals(model_semantic))
par(mfrow = c(1, 1))
dev.off()

# Figure ------------------------------------------------------------------
cell_plot <- as.data.frame(summary(cell_emm, infer = c(TRUE, TRUE))) %>%
  mutate(
    Tool_label = factor(Tool, levels = c("CAD", "PP"),
                        labels = c("CAD", "Pen-and-paper")),
    Task = factor(Task, levels = c("Bed", "Cabinet"))
  )

raw_plot <- df %>%
  mutate(Tool_label = factor(Tool, levels = c("CAD", "PP"),
                             labels = c("CAD", "Pen-and-paper")))

condition_colors <- c("CAD" = "#4C78A8", "Pen-and-paper" = "#F58518")
condition_shapes <- c("CAD" = 21, "Pen-and-paper" = 24)
condition_linetypes <- c("CAD" = "solid", "Pen-and-paper" = "dashed")
dodge <- position_dodge(width = 0.16)

semantic_plot <- ggplot() +
  geom_point(
    data = raw_plot,
    aes(Task, Semantic, shape = Tool_label, fill = Tool_label, color = Tool_label),
    position = position_jitterdodge(jitter.width = 0.06, dodge.width = 0.16),
    size = 2.2, alpha = 0.22
  ) +
  geom_line(
    data = cell_plot,
    aes(Task, emmean, color = Tool_label, linetype = Tool_label, group = Tool_label),
    position = dodge, linewidth = 0.95
  ) +
  geom_errorbar(
    data = cell_plot,
    aes(Task, ymin = lower.CL, ymax = upper.CL, color = Tool_label, group = Tool_label),
    position = dodge, width = 0.055, linewidth = 0.8
  ) +
  geom_point(
    data = cell_plot,
    aes(Task, emmean, shape = Tool_label, fill = Tool_label, color = Tool_label),
    position = dodge, size = 4.4
  ) +
  scale_color_manual(values = condition_colors, name = "Design medium") +
  scale_fill_manual(values = condition_colors, name = "Design medium") +
  scale_shape_manual(values = condition_shapes, name = "Design medium") +
  scale_linetype_manual(values = condition_linetypes, name = "Design medium") +
  labs(x = "Furniture task", y = "Semantic distance") +
  theme_classic(base_size = 15) +
  theme(legend.position = "top")

print(semantic_plot)

ggsave(file.path(figures_dir, "Semantic_Distance_Tool_Task.png"), semantic_plot,
       width = 7.5, height = 5.5, dpi = 600, bg = "white")
ggsave(file.path(figures_dir, "Semantic_Distance_Tool_Task.pdf"), semantic_plot,
       width = 7.5, height = 5.5)
