# Space Ranger Tutorial

## Installing 10x Space Ranger and Coda Environment

### 1. Installing 10x Space Ranger and Coda Environment

To install 10x Space Ranger and set up the corresponding Coda environment:

1. **Download Space Ranger:**
   - Go to the 10x Genomics website and download the latest version of Space Ranger: https://support.10xgenomics.com/spatial-gene-expression/software/downloads/latest
   - Follow the installation instructions specific to your operating system.

2. **Install Coda Environment:**
   - Create a conda environment for running Space Ranger:
     ```bash
     conda create -n spaceranger_env python=3.8
     conda activate spaceranger_env
     ```
   - Install dependencies for the Coda environment:
   - 
     ```yaml
      name: space-ranger
      channels:
        - conda-forge
        - bioconda
        - defaults
      dependencies:
        - python=3.7
        - numpy
        - pandas
        - matplotlib
        - scipy
        - scikit-learn
        - h5py
        - pyyaml
        - biopython
        - seaborn
        - zlib
        - libpng
        - curl
        - bcftools
        - samtools
        - bedtools
        - openjdk
        - scanpy
        - umap-learn
        - torch
     ```
    

3. **Verify Installation:**
   - To verify if Space Ranger was installed properly, run:
     ```bash
     spaceranger --version
     ```
---

### 2. Processing Files for Space Ranger

#### a. Convert GFF to GTF for Space Ranger

- Use `gffread` to convert the gene annotation file from GFF to GTF format:

```bash
gffread Gbi_Genes.gff -T -o Gbi_Genes.gtf
```

#### b. Count Contig Length

* To count the length of the contigs in the genome file , run the following command:

```bash
awk '/^>/ {if (seqlen > 0) print name, seqlen; name = substr($0, 2); seqlen = 0} 
     /^[^>]/ {seqlen += length($0)} 
     END {print name, seqlen}' file.fasta > contig_lengths.txt
```

#### c. Eliminate Lines Greater Than Contig Length

* Use the following command to correct any lines where the `end` value is greater than the contig length, and then update the `Gbi_Genes.gtf` file:

```bash
awk 'BEGIN {while ((getline < "contig_lengths.txt") > 0) {contig[$1] = $2}} 
     {seqname = $1; start = $4; end = $5; if (contig[seqname] && end > contig[seqname]) { $5 = contig[seqname] }} 
     {print $0}' Gbi_Genes.gtf > Gbi_Genes_clean.gtf && mv Gbi_Genes_clean.gtf Gbi_Genes.gtf
```

#### d. Eliminate Lines with Unmatched Columns

* Remove any lines with unmatched columns (empty values) from `Gbi_Genes.gtf` :

```bash
awk 'BEGIN {FS=OFS="\t"} {if (NF == 9) print $0}' Gbi_Genes.gtf > Gbi_Genes_clean.gtf && mv Gbi_Genes_clean.gtf Gbi_Genes.gtf
```

## Execution of Space Ranger

#### mkref

```bash
./spaceranger mkref \
  --genome=GENOME \
  --fasta=../GENOME PATH/Gbi_genome.fa \
  --genes=../GTF PATH/Gbi_Genes_clean.gtf
```

#### count

```bash
./spaceranger count \
  --id=VISIUM ID \
  --transcriptome=/TRANSCRIPTOME PATH/ \
  --image=../IMAGE PATH/image.tif \
  --fastqs=../FASTQ PATH/ \
  --unknown-slide=visium-1 \ 
  --create-bam=false
```

#### count multiple slides

```bash
# Sample A/B/C/D
~/projects/spaceranger-3.1.2/spaceranger count \
  --id=Sample_IAB \
  --transcriptome= transcriptome_path #/gbim1/Annotation/Genes/ # Repo output mkref (fasta genes star reference.json)
  --fastqs=/home/martinpl/projects/datashare/grillus/2/A \ # Repo with 2 FASTQs / sample
  --sample=GMRS238 \   #38/39/40/41
  --image= Image_path
  --slide=V13B16-384 \ #same for each
  --area=A1 \ #B1/C1/D1
  --create-bam=false

```

## Output overview

| Output                        | Description                                                                 |
|-------------------------------|-----------------------------------------------------------------------------|
| Analysis  | Results from clustering, differential expression analysis, and dimensionality reduction techniques (PCA, t-SNE, UMAP). |
| Deconvolution      | Results of deconvolution analysis using different values of k (from k=2 to k=8). |
| Filtered feature barcode matrix  | Processed gene expression matrix, where `matrix.mtx` contains gene expression data, `features.tsv` lists gene identifiers, and `barcodes.tsv` contains cell barcode identifiers. |
| Raw feature barcode matrix     | Unprocessed gene expression matrix containing raw data for features and barcodes. |
| Spatial                       | Spatial information about tissue spots, including spatial coordinates and gene expression data. |
| metrics_summary.csv            | A CSV file summarizing various metrics and statistics from the analysis process. |
| molecule_info.h5               | HDF5 file containing information about the molecules, such as read counts and molecule IDs. |
| web_summary.html               | HTML summary report of the analysis, including interactive visualizations and overall data insights. |





















