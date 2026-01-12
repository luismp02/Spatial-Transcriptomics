#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import numpy as np
import scanpy as sc
import warnings

warnings.filterwarnings("ignore")

# ==========================================
# Thresholds per sample (from your table)
# Format: sample_id: (min_count, max_count)
# ==========================================
thresholds = {
    "V43T08-041_A":  (900, np.inf),
    "V43T08-041_D":  (0, 50000),
    "V43T08-051_A":  (0, 90000),
    "V43T08-051_D":  (2500, np.inf),
    "V44L01-325_A":  (0, np.inf),
    "V44L01-325_D":  (0, np.inf),
    "V44L22-378_A":  (0, np.inf),
    "V44L22-378_D":  (0, np.inf),
    "V44L23-362_A":  (0, 39000),
    "V44L23-362_D":  (0, np.inf),
    "V44L23-391_A":  (0, np.inf),
    "V44L23-391_D":  (0, np.inf),
}

# Base paths
base_dir = "/home/martinpl/projects/datashare/PDAC/visium_PDAC"

slides = [
    "Visium_FFPE_V43T08-051_D", "Visium_FFPE_V44L22-378_D", "Visium_FFPE_V44L23-391_D",
    "Visium_FFPE_V43T08-041_A", "Visium_FFPE_V44L01-325_A", "Visium_FFPE_V44L23-362_A",
    "Visium_FFPE_V43T08-041_D", "Visium_FFPE_V44L01-325_D", "Visium_FFPE_V44L23-362_D",
    "Visium_FFPE_V43T08-051_A", "Visium_FFPE_V44L22-378_A", "Visium_FFPE_V44L23-391_A"
]

base_dir = "/home/martinpl/projects/datashare/PDAC/visium_PDAC"

# ==========================================
# Main loop
# ==========================================
for slide in slides:
    # extract short ID (after last underscore)
    short_id = slide.replace("Visium_FFPE_", "")
    if short_id not in thresholds:
        print(f"⚠️  No thresholds defined for {slide} (ID: {short_id}), skipping...")
        continue

    min_thr, max_thr = thresholds[short_id]
    print(f"\n--- Processing {slide} ({short_id}) with thresholds {min_thr} - {max_thr}")

    sample_input_dir = os.path.join(base_dir, slide, "outs")
    sample_output_dir = os.path.join(base_dir, slide)
    os.makedirs(sample_output_dir, exist_ok=True)

    if not os.path.exists(sample_input_dir):
        print(f"⚠️  Skipping {slide}, directory not found: {sample_input_dir}")
        continue

    # Load Visium data
    adata = sc.read_visium(sample_input_dir)

    # Compute QC metrics and filter spots by total counts
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    keep = (adata.obs["total_counts"] >= min_thr) & (adata.obs["total_counts"] <= max_thr)
    adata = adata[keep, :]

    # Normalize and log-transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Export to CSV (spots as rows, genes as columns)
    df = pd.DataFrame(
        adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X,
        index=adata.obs_names,
        columns=adata.var_names
    )

    csv_path = os.path.join(sample_output_dir, f"{slide}_filtered.csv")
    df.to_csv(csv_path)

    print(f"✓ Saved: {csv_path}")

