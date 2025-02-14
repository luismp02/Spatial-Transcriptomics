# PDAC datashare overview

## Data Directory
```
user@f-dahu:/bettik/PROJECTS/pr-epimed/martinpl/projects/datashare/
```

### Visium data
Primary datasets stored in the following path (bettik path to be updated):

```
user@f-dahu:/bettik/PROJECTS/pr-epimed/martinpl/projects/datashare/visium_PDAC
```

This directory contains 12 samples, each corresponding to a different experiment:

| Sample Directory | Description |
|------------------|-------------|
| Visium_FFPE_V43T08-051_D | FFPE sample ID V43T08-051 (D) |
| Visium_FFPE_V44L22-378_D | FFPE sample ID V44L22-378 (D) |
| Visium_FFPE_V44L23-391_D | FFPE sample ID V44L23-391 (D) |
| Visium_FFPE_V43T08-041_A | FFPE sample ID V43T08-041 (A) |
| Visium_FFPE_V44L01-325_A | FFPE sample ID V44L01-325 (A) |
| Visium_FFPE_V44L23-362_A | FFPE sample ID V44L23-362 (A) |
| Visium_FFPE_V43T08-041_D | FFPE sample ID V43T08-041 (D) |
| Visium_FFPE_V44L01-325_D | FFPE sample ID V44L01-325 (D) |
| Visium_FFPE_V44L23-362_D | FFPE sample ID V44L23-362 (D) |
| Visium_FFPE_V43T08-051_A | FFPE sample ID V43T08-051 (A) |
| Visium_FFPE_V44L22-378_A | FFPE sample ID V44L22-378 (A) |
| Visium_FFPE_V44L23-391_A | FFPE sample ID V44L23-391 (A) |

### Sample Structure

Each sample directory contains key output files organized in subdirectories:

### Filtered Feature Barcode Matrix
Path: `/outs/filtered_feature_bc_matrix/`

| File | Description |
|------|-------------|
| `barcodes.tsv.gz` | List of barcode identifiers |
| `features.tsv.gz` | List of detected features (genes) |
| `matrix.mtx.gz` | Matrix of feature counts per barcode |

### Spatial Data
Path: `/outs/spatial/`

| File | Description |
|------|-------------|
| `tissue_hires_image.png` | High-resolution image of the tissue sample |
| `scalefactors_json.json` | JSON file containing scaling factors for spatial mapping |
| `tissue_positions.csv` | CSV file mapping barcodes to spatial coordinates |


## Reference Data
* Reference data
```
user@f-dahu:/bettik/PROJECTS/pr-epimed/martinpl/projects/datashare/sc_PDAC
```
This directory contains two main datasets:

### GSE194247_RAW

| Sample ID | Files |
|-----------|-------|
| GSM5831620_5_GEX_4 | `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz` |
| GSM5831621_5_GEX_5 | `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz` |
| GSM5831622_5_GEX_6 | `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz` |
| GSM5831623_5_GEX_9 | `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz` |
| GSM5831624_GEX_45_MM | `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz` |

### GSE235449_RAW

| Sample ID | Files |
|-----------|-------|
| GSM7502530 | `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz` |


