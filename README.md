# EEG-study-of-cognitive-load-and-design-fixation-in-CAD-versus-pen-and-paper-conceptual-design.

This repository contains the preprocessing, feature-extraction, statistical-analysis code, and de-identified derived data used to support the manuscript:

**Does design medium affect cognitive load and design fixation during divergent thinking? A mixed-methods study**

The study compares conceptual design using **computer-aided design (CAD)** and **pen-and-paper (PP)** in a within-subject experiment. Behavioral, self-report, and electroencephalography (EEG) measures were used to examine two primary outcomes:

1. **Cognitive load**

   * NASA Task Load Index (NASA-TLX)
   * EEG-based Cognitive Load Index (CLI)

2. **Design fixation**

   * Semantic distance of generated ideas
   * EEG task-related power (TRP), ERD/ERS patterns, and multiple correspondence analysis (MCA)

---

## Repository Structure

```text
CHB_reproducibility_code_final/
│
├── analysis/
│   ├── 03_NASA_TLX_analysis.R
│   ├── 04_CLI_analysis.R
│   ├── 05_Semantic_Distance_analysis.R
│   ├── 06_TRP_preparation.R
│   ├── 07_MCA_analysis.R
│   └── 08_EEG_quality_summary.R
│
├── data/
│   ├── PSD_CLI.csv
│   ├── MCA_ERS_ERD.csv
│   ├── NASA_semantic.csv
│   └── README.md
│
├── preprocessing/
│   ├── 00_Segment_EEG.py
│   ├── 01_Preprocessing_reference_global.m
│   └── 02_Feature_Extraction_as_run.m
│
└── README.md
```

---

## Data

All shared datasets are de-identified derived datasets used for the analyses reported in the manuscript.

Participant identifiers are standardized consistently across datasets as `PTP01`–`PTP28`.

### `PSD_CLI.csv`

Contains participant-level EEG features organized by participant, design condition, experimental phase, and furniture task.

The dataset contains derived EEG measures generated from the feature-extraction pipeline, including:

* power spectral density (PSD) measures
* frequency-band power measures
* EEG-based Cognitive Load Index (CLI)

This dataset is used primarily by:

```text
04_CLI_analysis.R
```

and provides EEG-derived measures used in subsequent analyses.

---

### `NASA_semantic.csv`

Contains the de-identified behavioral and questionnaire-derived measures used in the cognitive-load and behavioral design-fixation analyses.

The dataset includes information used for:

* NASA-TLX workload analysis
* semantic-distance analysis
* design condition
* furniture task
* participant-level experimental information required by the corresponding models

This dataset is used primarily by:

```text
03_NASA_TLX_analysis.R
05_Semantic_Distance_analysis.R
```

---

### `MCA_ERS_ERD.csv`

Contains the categorized EEG responses used for the ERD/ERS and multiple correspondence analysis.

EEG features are categorized into:

* ERD
* No Change
* ERS

These categorical EEG states are used to examine multivariate neural patterns associated with CAD and pen-and-paper design.

This dataset is used by:

```text
07_MCA_analysis.R
```

---

### Raw EEG Data

Raw participant-level EEG recordings are **not included in this repository** because public sharing of the raw recordings was not covered by participant consent and the ethical approval for the study.

The repository instead provides de-identified derived EEG measures and analysis code needed to reproduce the statistical analyses reported in the manuscript.

---

# EEG Processing Pipeline

The EEG processing workflow is contained in the `preprocessing/` folder.

## 1. EEG Segmentation

```text
preprocessing/00_Segment_EEG.py
```

Segments the continuous EEG recordings according to participant, design condition, and experimental phase.

The resulting segmented files provide the input to the EEG preprocessing pipeline.

---

## 2. EEG Preprocessing

```text
preprocessing/01_Preprocessing_reference_global.m
```

EEG preprocessing was performed in MATLAB using EEGLAB.

The major processing steps include:

* 0.5–40 Hz band-pass filtering
* line-noise filtering
* automated bad-channel and bad-segment detection using `clean_rawdata`
* Artifact Subspace Reconstruction (ASR)
* interpolation of removed channels
* common-average re-referencing
* Independent Component Analysis using Picard
* automated component classification using ICLabel
* removal of artifactual independent components

The preprocessing script represents the pipeline applied to the EEG data used in the manuscript.

---

## 3. EEG Feature Extraction

```text
preprocessing/02_Feature_Extraction_as_run.m
```

Extracts frequency-domain EEG features from the cleaned EEG recordings.

The script calculates:

* power spectral density
* frequency-band power
* frontal theta power
* parietal alpha power
* EEG-based Cognitive Load Index (CLI)

The resulting participant-level EEG measures are included in:

```text
data/PSD_CLI.csv
```

---

# Statistical Analyses

Analysis scripts are located in the `analysis/` folder.

## `03_NASA_TLX_analysis.R`

Analyzes subjective workload differences between CAD and pen-and-paper using data from:

```text
data/NASA_semantic.csv
```

The script includes:

* descriptive statistics
* paired-sample comparisons for the six NASA-TLX dimensions
* Bonferroni correction
* paired-sample effect sizes
* manuscript figure generation

The six NASA-TLX dimensions are:

* Mental Demand
* Physical Demand
* Temporal Demand
* Performance
* Effort
* Frustration

---

## `04_CLI_analysis.R`

Analyzes the EEG-based Cognitive Load Index using:

```text
data/PSD_CLI.csv
```

Baseline correction is calculated as:

```text
Delta CLI = Free-design CLI - Baseline CLI
```

The script includes:

* creation of complete CAD–PP participant pairs
* descriptive statistics
* assessment of paired difference scores
* paired-sample t-test
* Cohen's dz
* Wilcoxon signed-rank sensitivity analysis
* manuscript figure generation

---

## `05_Semantic_Distance_analysis.R`

Analyzes behavioral design fixation using semantic-distance measures from:

```text
data/NASA_semantic.csv
```

Higher semantic distance represents greater divergence from the original furniture prompt, whereas lower semantic distance represents greater semantic similarity to the original object.

The primary model evaluates:

```text
Semantic Distance ~ Design Tool × Furniture Task + (1 | Participant)
```

where:

* Design Tool = CAD vs. pen-and-paper
* Furniture Task = bed vs. cabinet
* Participant = random intercept

The script includes:

* linear mixed-effects modeling
* estimated marginal means
* CAD vs. pen-and-paper contrast
* furniture-task contrast
* effect-size estimation
* basic model diagnostics
* manuscript figure generation

---

## `06_TRP_preparation.R`

Calculates baseline-corrected task-related power (TRP) measures for the EEG design-fixation analyses.

TRP is calculated from log-transformed spectral power as:

```text
TRP = log(Task Power) - log(Baseline Power)
```

The resulting frequency- and region-specific EEG features are used in the subsequent ERD/ERS and MCA analyses.

---

## `07_MCA_analysis.R`

Performs the multivariate EEG analysis of fixation-related neural activity.

The workflow includes:

* categorization of baseline-corrected TRP values into ERD, No Change, and ERS
* Multiple Correspondence Analysis (MCA)
* extraction of MCA dimensions
* category coordinates
* contributions to MCA dimensions
* supplementary feature-level comparisons
* manuscript figure generation

---

## `08_EEG_quality_summary.R`

Summarizes EEG data retention after preprocessing.

The script calculates retained-data proportions by experimental condition and supports the EEG quality-control information reported in the manuscript.

---

# Recommended Analysis Order

```text
03_NASA_TLX_analysis.R
        ↓
04_CLI_analysis.R
        ↓
05_Semantic_Distance_analysis.R

06_TRP_preparation.R
        ↓
07_MCA_analysis.R

08_EEG_quality_summary.R
```

The scripts are numbered for organization rather than because all analyses depend sequentially on one another.

---

# Software

The analyses were conducted using Python, MATLAB/EEGLAB, and R.

## Python

Used for EEG segmentation.

```text
Python 3.x
```

## MATLAB / EEGLAB

Used for EEG preprocessing and feature extraction.

Major EEGLAB functions and plugins used include:

* EEGLAB
* `clean_rawdata`
* Artifact Subspace Reconstruction (ASR)
* Picard ICA
* ICLabel

## R

Used for statistical analyses and visualization.

Major packages include:

```text
tidyverse
effectsize
lmerTest
emmeans
afex
FactoMineR
factoextra
ggplot2
```

Individual scripts specify the packages required for each analysis.

---

# Reproducibility

This repository provides the code and de-identified derived data required to reproduce the primary statistical analyses reported in the associated manuscript.

Because raw human-participant EEG recordings cannot be publicly distributed under the existing participant consent and ethical approval, reproduction of the reported statistical results begins from the provided derived datasets.

Participant identifiers have been standardized as:

```text
PTP01–PTP28
```

No directly identifying participant information is included in the shared datasets.

---

# Data Availability

Raw participant-level EEG recordings and other potentially identifiable source data are not publicly available because public sharing of these data was not covered by participant consent or the ethical approval for this study.

De-identified derived data supporting the analyses reported in the manuscript, together with the corresponding analysis code, are provided in this repository.

---

# Citation

If using this repository, please cite the associated manuscript.

Citation information will be added following publication.

---

# Contact

Questions regarding the code or analyses may be submitted through the repository's GitHub Issues page.
