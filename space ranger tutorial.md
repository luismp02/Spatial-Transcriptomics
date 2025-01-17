# Space Ranger Pipeline Tutorial

This tutorial covers key steps for using the Space Ranger pipeline. Each section provides code examples and explanations to guide you through the process. 

---

## Table of Contents
1. [Creating Reference file with mkref](#mkref)
2. [Counting with count](#count)
3. [Outputs](#outputs)
---
## Creating Reference file with mkref

#### Transform GTT to GTF
GTF format compatible with Space Ranger:

```bash
gffread Genes.gff -T -o Genes.gtf
```

---

#### Count Contig Length

Calculate the length of each contig in the genome FASTA file:

```bash
awk '/^>/ {if (seqlen > 0) print name, seqlen; name = substr($0, 2); seqlen = 0}
    /^[^>]/ {seqlen += length($0)}
    END {print name, seqlen}' file.fasta > contig_lengths.txt
```

---

##### Clean GTF File

Adjust entries in the GTF file where the end position exceeds the contig length:

```bash
awk 'BEGIN {
    while ((getline < "contig_lengths.txt") > 0) {
        contig[$1] = $2
    }
}
{
    seqname = $1
    start = $4
    end = $5
    if (contig[seqname] && end > contig[seqname]) {
        $5 = contig[seqname]
    }
    print $0
}' Genes.gtf > Genes_clean.gtf && mv Genes_clean.gtf Genes.gtf
```

Ensure all rows in the GTF file have the correct number of columns (9 in this case):

```bash
awk 'BEGIN {FS=OFS="\t"} {if (NF == 9) print $0}' Genes.gtf > Genes_clean.gtf && mv Genes_clean.gtf Genes.gtf
```

---

Run the `mkref` command to generate a reference for Space Ranger:

```bash
./spaceranger mkref \
  --genome=GRILLUS_GENOME \
  --fasta=../path/genome.fa \
  --genes=../path/Genes_clean.gtf
```

## Counting with Space Ranger
Input: FASTQ folder path, slide ID or slide model, image.tif
Run the `count` function to align reads and generate spatial expression data. Replace `slide id` with the appropriate identifier if available:

```bash
./spaceranger count \
  --id=GRILLUS_ST7 \
  --transcriptome=/path/spaceranger-3.1.1/GENOME/ \
  --image=../path/image.tif \
  --fastqs=../path/ \
  --unknown-slide=visium-1 \
  --create-bam=false
```

---

## Outputs 

The following outputs will be generated during the Space Ranger analysis:

1. **Filtered Feature Barcode Matrix**
   - `matrix.mtx`
   - `features.tsv`
   - `barcodes.tsv`
2. **Raw Feature Barcode Matrix**
3. **Spatial Data**
   - Processed image files and spatial coordinates.
   - Detected tissue image in jpeg format
   - Hires and Lowres image in png format
4. **Aditional files**
   - `metrics_summary.csv`
   - `molecule_info.h5`
   - `web_summary.html`

5. **Additional analyses**
- clustering
- differential expression
- PCA
- t-SNE
- UMAP
- deconvolution 

---

This tutorial provides a structured approach to processing spatial transcriptomics
