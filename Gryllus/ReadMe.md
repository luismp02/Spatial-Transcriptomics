# Gryllus – Spatial Transcriptomics Workflows

This directory contains workflows and analysis scripts related to Gryllus bimaculatus spatial transcriptomics.

Unlike PDAC, Gryllus analyses are more exploratory and organism-specific, and the pipeline reflects this distinction.

## Overview

The Gryllus section is organized around two complementary layers:

1. Infrastructure and preprocessing (Nextflow)  
2. Exploratory and biological analyses (standalone scripts)  

Only the infrastructure layer is formalized using Nextflow.

## Directory Structure

Gryllus/  
├── scripts/  
│   ├── load_data_nf.py  
│   ├── load_data_nf.R  
│   ├── spaceranger_mkref.sh  
│   ├── spaceranger_count.sh  
│   ├── 01_load_data.py  
│   ├── 02_gene_analysis.py  
│   ├── chrX-auto.py  
│   └── ucell.Rmd  
│  
└── nextflow/  
    ├── main.nf  
    ├── nextflow.config  
    ├── modules/  
    │   ├── spaceranger_mkref/  
    │   ├── spaceranger_count/  
    │   ├── load_py/  
    │   └── load_r/  
    └── assets/  
        └── samplesheet.tsv  

## Space Ranger Workflows

Location: Gryllus/scripts/ and Gryllus/nextflow/

### mkref

Builds a Space Ranger transcriptome reference using spaceranger mkref.  
The process is fully parameterized and can be run either standalone or through Nextflow.

### count

Runs spaceranger count for one or multiple Visium slides.  
Execution is driven by a samplesheet containing FASTQs, image paths, slide IDs, and capture areas.

## Data Loading

### Python loader

Script: load_data_nf.py

Converts Space Ranger outputs into an AnnData (.h5ad) object, annotating genes with regeneration, growth, and scarring signatures while preserving spatial metadata.

### R loader

Script: load_data_nf.R

Converts Space Ranger outputs into a Seurat (.rds) object, mirroring the Python logic for cross-language consistency.

Both loaders are integrated into Nextflow but remain fully usable as standalone scripts.

## Exploratory Analysis Scripts

Location: Gryllus/scripts/

Scripts such as 01_load_data.py, 02_gene_analysis.py, chrX-auto.py, and ucell.Rmd are not part of the Nextflow pipeline.

They represent exploratory analyses, hypothesis-driven investigations, and organism-specific biological questions.  
These scripts are intentionally kept outside the formal workflow.

## Nextflow Pipeline

Location: Gryllus/nextflow/

Supported modes include mkref, count, load_py, load_r, load_only, and all.

Example usage:

nextflow run main.nf --mode count

The pipeline is designed to standardize heavy preprocessing steps while remaining flexible for downstream exploratory analysis.

## Design Philosophy

Only repetitive, infrastructure-heavy steps are automated.  
Biological interpretation remains script-driven.  
Outputs and environments are excluded from version control.  
The pipeline prioritizes clarity over maximal automation.
