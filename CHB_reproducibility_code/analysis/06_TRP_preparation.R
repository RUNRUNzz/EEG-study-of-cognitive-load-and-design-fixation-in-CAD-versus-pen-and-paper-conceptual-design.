# Task-related power (TRP) preparation for the CHB fixation analysis
# Computes free-design minus modality-matched baseline log power.

library(tidyverse)

# Paths -------------------------------------------------------------------
data_file <- file.path("data", "Frequency_bands_power.csv")
results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE)

bands <- read_csv(data_file, show_col_types = FALSE) %>%
  mutate(
    Participant = trimws(as.character(Participant)),
    Condition = toupper(trimws(as.character(Condition))),
    Phase = tolower(trimws(as.character(Phase))),
    Channel = toupper(trimws(as.character(Channel)))
  ) %>%
  filter(Condition %in% c("CAD", "PP"), Phase %in% c("baseline", "free_design"))

# Channel regions used in the original TRP analysis -----------------------
regions <- list(
  Entire = c("FP1", "F7", "F8", "T4", "T6", "T5", "T3", "FP2", "O1", "P3",
             "PZ", "F3", "FZ", "F4", "C4", "P4", "POZ", "C3", "CZ", "O2"),
  Frontal = c("FP1", "FP2", "F7", "F3", "FZ", "F4", "F8"),
  Parietal = c("P3", "PZ", "P4"),
  Occipital = c("O1", "O2", "POZ"),
  FZ = "FZ",
  PZ = "PZ"
)

region_map <- enframe(regions, name = "Region", value = "Channel") %>%
  unnest(Channel)

# Convert band powers to long format and log-transform --------------------
long_power <- bands %>%
  pivot_longer(Delta:Gamma, names_to = "Band", values_to = "Power") %>%
  filter(is.finite(Power), Power > 0) %>%
  mutate(logPower = log(Power))

# Match free-design power to the corresponding baseline ------------------
baseline <- long_power %>%
  filter(Phase == "baseline") %>%
  select(Participant, Condition, Channel, Band, baseline_logPower = logPower)

free_design <- long_power %>%
  filter(Phase == "free_design") %>%
  select(Participant, Condition, Channel, Band, task_logPower = logPower)

trp_channel <- free_design %>%
  inner_join(baseline, by = c("Participant", "Condition", "Channel", "Band")) %>%
  mutate(TRP = task_logPower - baseline_logPower)

# Average channels within each region -------------------------------------
trp_region <- trp_channel %>%
  inner_join(region_map, by = "Channel") %>%
  group_by(Participant, Condition, Region, Band) %>%
  summarise(TRP = mean(TRP, na.rm = TRUE), .groups = "drop")

# Features used in the original MCA input --------------------------------
trp_features <- trp_region %>%
  filter(
    (Band == "Theta" & Region %in% c("Entire", "Frontal")) |
      (Band == "Alpha" & Region %in% c("Entire", "Parietal", "Occipital", "PZ")) |
      (Band == "Beta" & Region == "Entire")
  ) %>%
  transmute(
    Participant,
    Condition,
    Metric = paste0("TRP_", Band, "_", Region, "_diff_FD"),
    Value = TRP
  )

# Theta/alpha TRP ratio for the whole scalp -------------------------------
theta_alpha <- trp_region %>%
  filter(Region == "Entire", Band %in% c("Theta", "Alpha")) %>%
  select(Participant, Condition, Band, TRP) %>%
  pivot_wider(names_from = Band, values_from = TRP) %>%
  transmute(
    Participant,
    Condition,
    Metric = "TRP_log(Theta/Alpha)_diff_FD",
    Value = Theta - Alpha
  )

# Fz-theta / Pz-alpha TRP ratio -------------------------------------------
fz_pz <- trp_region %>%
  filter((Region == "FZ" & Band == "Theta") |
           (Region == "PZ" & Band == "Alpha")) %>%
  mutate(key = paste0(Region, "_", Band)) %>%
  select(Participant, Condition, key, TRP) %>%
  pivot_wider(names_from = key, values_from = TRP) %>%
  drop_na(FZ_Theta, PZ_Alpha) %>%
  transmute(
    Participant,
    Condition,
    Metric = "TRP_log(Theta_Fz/Alpha_Pz)_diff_FD",
    Value = FZ_Theta - PZ_Alpha
  )

trp_long <- bind_rows(trp_features, theta_alpha, fz_pz) %>%
  arrange(Participant, Condition, Metric)

trp_wide <- trp_long %>%
  pivot_wider(names_from = Metric, values_from = Value)

write_csv(trp_long, file.path(results_dir, "EEG_TRP_free_design_long.csv"))
write_csv(trp_wide, file.path(results_dir, "EEG_TRP_free_design_wide.csv"))

print(head(trp_wide))
