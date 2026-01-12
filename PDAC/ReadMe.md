# PDAC – Spatial Transcriptomics Pipeline

This directory contains the PDAC (Pancreatic Ductal Adenocarcinoma) spatial transcriptomics workflows developed in this project.  
The focus is on reproducible preprocessing, deconvolution, and simulation, implemented using Nextflow while keeping all scientific scripts executable outside of the pipeline.

## Overview

The PDAC analysis is organized into three largely independent components:

1. Preprocessing  
2. Deconvolution (SPOTlight & SPADE)  
3. Simulation (scCube)  

Additionally, clustering is handled in a separate, external pipeline and is intentionally not coupled to the main PDAC workflow.

## Directory Structure

PDAC/  
├── Preprocessing/  
│   └── scripts/  
│       └── preprocess_integrated.R  
│  
├── deconvolution/  
│   ├── SPOTlight/  
│   │   └── scripts/  
│   │       └── SPOTlight.R  
│   │  
│   └── SPADE/  
│       ├── scripts/  
│       │   └── SPADE.R  
│       └── SpaGCN/  
│           └── SpaGCN and domain clustering code  
│  
├── simulations/  
│   └── simulation_scCube/  
│       └── scCube scripts only (no outputs or environments)  
│  
└── nextflow/  
    ├── main.nf  
    ├── nextflow.config  
    ├── modules/  
    │   ├── preprocess/  
    │   ├── spotlight/  
    │   ├── spade/  
    │   └── simulation/  
    └── assets/  
        └── slides.txt  

## Preprocessing

Location: PDAC/Preprocessing/scripts/

The preprocessing step prepares the scRNA-seq reference and the Visium spatial data.  
It performs normalization, dimensionality reduction, and annotation transfer.

This step is always executed first in the Nextflow pipeline and serves as a logical dependency for downstream analyses.

## Deconvolution

### SPOTlight

Location: PDAC/deconvolution/SPOTlight/scripts/

SPOTlight performs topic-based deconvolution of Visium spots using a scRNA-seq reference with fine-grained annotations.  
It is implemented as a single-slide script (SPOTlight.R) with command-line arguments and integrated into Nextflow via a dedicated module.

### SPADE

Location: PDAC/deconvolution/SPADE/

SPADE uses spatial domain information (e.g. SpaGCN output) to perform domain-aware deconvolution.  
The workflow includes fused reference construction and domain-specific modeling.

SPOTlight and SPADE are parallel alternatives, not sequential steps.

## Simulation (scCube)

Location: PDAC/simulations/simulation_scCube/

Simulation is intentionally independent from preprocessing and deconvolution.  
It is based on scCube and is mainly used for pedagogical exploration and method testing.

Simulation is integrated into Nextflow as an optional branch.  
All outputs and environments are excluded from version control.

## Clustering (Important Note)

Clustering is not part of this pipeline.

The clustering workflow (GraphST / SpaGCN-based) has its own Nextflow pipeline, its own documentation, and its own repository structure.  
PDAC deconvolution workflows consume clustering results but do not generate them.

## Nextflow Pipeline

Location: PDAC/nextflow/

The pipeline supports different execution modes, including preprocess, spotlight, spade, simulation, and all.

Example usage:

nextflow run main.nf --mode spotlight

All outputs are written to a configurable results directory and are excluded from version control.

## Design Philosophy

Scripts remain executable outside Nextflow.  
Nextflow is used for orchestration rather than hiding logic.  
Outputs are never versioned.  
Large environments and binary data are excluded.  
Each methodological block remains conceptually independent.
