"""Segment raw EEG CSV files into baseline, free-design, and time-constrained phases.

Edit PROJECT_ROOT before running. The script expects the original study naming
convention for EEG and annotation CSV files.
"""

from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path("/path/to/EEG_2")
RAW_DIR = PROJECT_ROOT / "0_raw"
MARKER_DIR = PROJECT_ROOT / "0_marker"
OUTPUT_DIR = PROJECT_ROOT / "1_segmented"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

EEG_COLUMNS = [
    "Timestamp", "Fp1", "F7", "F8", "T4", "T6", "T5", "T3",
    "Fp2", "O1", "P3", "Pz", "F3", "Fz", "F4", "C4",
    "P4", "POz", "C3", "Cz", "O2",
]
EEG_CHANNELS = [c for c in EEG_COLUMNS if c != "Timestamp"]

for eeg_file in sorted(RAW_DIR.glob("*.csv")):
    print(f"Processing: {eeg_file.name}")

    parts = eeg_file.stem.split(" ")
    if len(parts) < 2:
        raise ValueError(f"Unexpected EEG filename: {eeg_file.name}")

    participant, condition = parts[0], parts[1]

    eeg = pd.read_csv(eeg_file, header=25)
    eeg = eeg[EEG_COLUMNS]
    eeg = eeg.dropna(subset=EEG_CHANNELS, how="all")

    marker_file = MARKER_DIR / f"{participant} {condition} Individual Annotation Data.csv"
    markers = pd.read_csv(marker_file, header=2)

    for _, row in markers.iterrows():
        phase = str(row["Marker Name"]).strip().lower().replace(" ", "_")
        start_ms = float(row["Start Time (ms)"])
        end_ms = float(row["End Time (ms)"])

        segment = eeg.loc[
            (eeg["Timestamp"] >= start_ms) & (eeg["Timestamp"] <= end_ms),
            EEG_CHANNELS,
        ]

        output_file = OUTPUT_DIR / f"{participant}_{condition}_{phase}.csv"
        segment.to_csv(output_file, index=False)
        print(f"  Saved: {output_file.name}")
