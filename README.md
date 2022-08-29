## MISET-seq
### The first high-throughput spatial co-omics method, named as microfluidic indexing based spatial epigenome and transcriptome sequencing (MISET-seq).

![image](https://raw.githubusercontent.com/gpenglab/MISET-seq/main/MISET-seq.png)

### MISET-seq data analysis pipeline
### 1. Download and install cellranger-arc software from 10x genomics and replace their default barcode file:
a. Enter cellranger-arc-2.0.1/lib/python/cellranger/barcodes/ and cellranger-arc-2.0.1/lib/python/atac/barcodes/ directory;

b. Replace default barcode file "737K-arc-v1.txt.gz" with your new custom barcode (in barcode file, also named as "737K-arc-v1.txt.gz").

### 2. Generate ATAC fragments file and gene expression matrix file from MISET-seq data:
bash MISET-seq.sh

### 3. After fill tissue region with white and non-tissue with black in photoshop or other image processing software manually:
python Grid_filter.py

### 4. Analysis MISET-seq data in ArchR:
Rscript MISET-seq-ArchR.R

### 5. Analysis MISET-seq data in Signac:
Rscript MISET-seq-Signac.R
