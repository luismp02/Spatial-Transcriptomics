# Spatial-Transcriptomics

This repository contains workflows, tutorials, and scripts for two projects:

- **PDAC**: Analysis of pancreatic ductal adenocarcinoma (PDAC) using spatial transcriptomics.
- **Gryllus**: Spatial transcriptomics applied to the study of tissue regeneration in cricket species.

Each project contains:

- **Tutorials**: for installing and using ST tools (10xSpaceRanger, cell2location, GPU cluster)
- **Scripts**: DONE: loading data, integrating data, processing images, and integrating single-cell reference datasets.
- **Data Information**: Details on datashare of each project

## Repository Structure

```
spatial-transcriptomics-project/
│-- README.md
│-- PDAC/
│   │-- README.md
│   │-- data_info.md
│   │-- tutorials/
│   │   │-- installing_cell2location.md
│   │-- scripts/
│   │   │-- 01_load_data.py
│   │   │-- 02_load_sc_reference.py
│   │   │-- 03_exploratory_analysis
│
│-- Gryllus/
│   │-- README.md
│   │-- data_info.md
│   │-- tutorials/
│   │   │-- space_ranger.md
│   │   │-- data_integration.md
│   │-- scripts/
│   │   │-- 01_load_data.py
│   │   │-- 02_integration.py
│   │   │-- 03_gene_info_analysis.py
```

**Clone the repository**:
   ```bash
   git clone https://github.com/yourusername/spatial-transcriptomics-project.git
   cd spatial-transcriptomics-project
   ```

## License


