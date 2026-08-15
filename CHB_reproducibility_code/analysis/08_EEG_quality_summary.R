# EEG retained-data summary used in the manuscript
# Input: one row per participant x condition x phase after preprocessing.

library(tidyverse)

data_file <- file.path("data", "preprocessing_quality_report.csv")
results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE)

quality <- read_csv(data_file, show_col_types = FALSE) %>%
  transmute(
    Participant = trimws(as.character(Participant)),
    Condition = toupper(trimws(as.character(Condition))),
    Phase = tolower(trimws(as.character(Phase))),
    Retained_Percent = as.numeric(Total_Data_Retained_Pct)
  ) %>%
  filter(Condition %in% c("CAD", "PP"))

# The manuscript values were summarized across the 84 phase recordings
# available in each condition (28 participants x 3 phases).
by_condition <- quality %>%
  group_by(Condition) %>%
  summarise(
    N = sum(!is.na(Retained_Percent)),
    Mean_Retention = mean(Retained_Percent, na.rm = TRUE),
    SD_Retention = sd(Retained_Percent, na.rm = TRUE),
    Min_Retention = min(Retained_Percent, na.rm = TRUE),
    Max_Retention = max(Retained_Percent, na.rm = TRUE),
    .groups = "drop"
  )

by_phase_condition <- quality %>%
  group_by(Condition, Phase) %>%
  summarise(
    N = sum(!is.na(Retained_Percent)),
    Mean_Retention = mean(Retained_Percent, na.rm = TRUE),
    SD_Retention = sd(Retained_Percent, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(by_condition, file.path(results_dir, "EEG_retention_by_condition.csv"))
write_csv(by_phase_condition, file.path(results_dir, "EEG_retention_by_phase_condition.csv"))

print(by_condition)
print(by_phase_condition)
