# Spatial-Transcriptomics

This repository gathers workflows, pipelines, tutorials, and analysis scripts developed around **spatial transcriptomics**, with a focus on two independent biological projects:

- **PDAC**: Spatial transcriptomics analysis of pancreatic ductal adenocarcinoma.
- **Gryllus**: Spatial transcriptomics applied to tissue regeneration in the cricket *Gryllus bimaculatus*.

The repository is designed as a **research and pedagogical resource**

## Project Structure

The repository is organized by biological project:

Spatial-Transcriptomics/  
├── PDAC/  
│   ├── Preprocessing/  
│   ├── deconvolution/  
│   ├── simulations/  
│   └── nextflow/  
│  
└── Gryllus/  
    ├── scripts/  
    └── nextflow/  

Each project contains its own documentation (`README.md`) describing its internal structure, workflows, and design choices.

## PDAC Project

The PDAC project focuses on spatial transcriptomics analysis of pancreatic ductal adenocarcinoma samples.

Main components include:

- Preprocessing of scRNA-seq references and Visium spatial data  
- Spatial deconvolution using **SPOTlight** and **SPADE**  
- Optional spatial data simulation using **scCube**  
- Workflow orchestration using **Nextflow**

Clustering is intentionally handled in a **separate, external pipeline** with its own documentation and is not coupled to the main PDAC workflow.

All heavy outputs, environments, and intermediate results are excluded from version control.

See `PDAC/README.md` for full details.

## Gryllus Project

The Gryllus project explores spatial transcriptomics in a non-model organism (*Gryllus bimaculatus*), with an emphasis on tissue regeneration.

It is structured around two layers:

- **Infrastructure and preprocessing**, automated using Nextflow  
- **Exploratory and biological analyses**, implemented as standalone scripts

Key features include:

- Space Ranger workflows (mkref and count)  
- Data loading into AnnData (Python) and Seurat (R) formats  
- Exploratory analyses on gene signatures, regeneration markers, and chromosome-level expression patterns

Only repetitive, infrastructure-heavy steps are formalized as pipelines; biological interpretation remains script-driven.

See `Gryllus/README.md` for full details.

## Getting Started

Clone the repository:

git clone https://github.com/luismp02/Spatial-Transcriptomics.git  
cd Spatial-Transcriptomics

Each project can then be explored independently by following the documentation in its respective directory.

## Status

This repository is under active development.  
Some scripts and workflows are marked as in progress and may evolve as the projects advance.

Contributions, reuse, and adaptation for research and teaching purposes are encouraged.
