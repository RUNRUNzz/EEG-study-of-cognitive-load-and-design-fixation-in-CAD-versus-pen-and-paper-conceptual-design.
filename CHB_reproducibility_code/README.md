# CHB reproducibility code

Cleaned and de-identified analysis code for the manuscript comparing CAD and pen-and-paper during divergent conceptual design.

## What is included

The package contains only code needed to support the analyses reported in the current manuscript:

1. EEG segmentation
2. EEG preprocessing
3. EEG band-power and CLI feature extraction
4. NASA-TLX paired comparisons
5. Baseline-corrected CLI comparison
6. Semantic-distance mixed-effects model
7. TRP preparation
8. MCA and feature-level follow-ups
9. EEG retained-data summary

Older exploratory analyses (FAA, connectivity, NASA-EEG correlations, FD/TC ANOVAs, and unrelated plots) were removed because they are not part of the current CHB results.

## Run order

Run from the project/repository root:

1. `preprocessing/00_Segment_EEG.py`
2. `preprocessing/01_Preprocessing_reference_global.m`
3. `preprocessing/02_Feature_Extraction_as_run.m`
4. `analysis/06_TRP_preparation.R`
5. `analysis/03_NASA_TLX_analysis.R`
6. `analysis/04_CLI_analysis.R`
7. `analysis/05_Semantic_Distance_analysis.R`
8. `analysis/07_MCA_analysis.R`
9. `analysis/08_EEG_quality_summary.R`

NASA-TLX and semantic-distance analyses can be run independently if their de-identified derived data are already available.

## De-identification

No author name or personal computer username is used in the cleaned analysis code. Generic `/path/to/...` placeholders are used for MATLAB/Python preprocessing paths.

Before public release, replace participant identifiers in all shared derived data with neutral IDs such as `S01` to `S28`. Keep the private ID mapping outside the public repository.

Do not upload raw EEG, audio, full verbal transcripts, consent forms, contact information, or other identifiable source data unless public sharing is explicitly covered by the study consent and ethics approval.

## Expected derived data files

Place these in `data/`:

- `nasa_tlx.csv`: `Participant`, `Condition`, `NASA_1` ... `NASA_6`
- `eeg_features.csv`: `Participant`, `Condition`, `Phase`, `CLI` (plus PSD columns if retained)
- `semantic_distance.csv`: `Participant`, `Tool`, `Task`, `Semantic`
- `Frequency_bands_power.csv`: `Participant`, `Condition`, `Phase`, `Channel`, `Delta`, `Theta`, `Alpha`, `Beta`, `Gamma`
- `preprocessing_quality_report.csv`: `Participant`, `Condition`, `Phase`, `Total_Data_Retained_Pct` (other QC columns may also be present)

## Software/packages

Python segmentation uses `pandas`.

MATLAB preprocessing/feature extraction requires EEGLAB and the plugins/functions used in the original pipeline, including clean_rawdata/ASR, Picard ICA, and ICLabel.

R analyses use: `tidyverse`, `effectsize`, `patchwork`, `lmerTest`, `emmeans`, `FactoMineR`, and `factoextra`.

## Important

Read `VALIDATION_NOTES.md` before making this repository public. The uploaded preprocessing file is the same version as the earlier preprocessing script and closely matches the manuscript preprocessing description, but one final rejection step is described incorrectly in the manuscript. In addition, the feature-extraction and MCA code contain settings that must be reconciled with the manuscript wording before final public release.
