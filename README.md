## MISAR-seq
### The first high-throughput spatial co-omics method, named as microfluidic indexing based spatial ATAC and RNA sequencing (MISAR-seq).

![image](https://raw.githubusercontent.com/gpenglab/MISAR-seq/main/MISAR-seq.png)

### MISAR-seq data analysis pipeline
### 1. Download and install cellranger-arc software from 10x genomics and replace their default barcode file:
a. Enter cellranger-arc-2.0.1/lib/python/cellranger/barcodes/ and cellranger-arc-2.0.1/lib/python/atac/barcodes/ directory;

b. Replace default barcode file "737K-arc-v1.txt.gz" with your new custom barcode (in Barcode, also named as "737K-arc-v1.txt.gz").

### 2. Generate ATAC fragments file and gene expression matrix file from MISAR-seq data:
bash MISAR-seq.sh

### 3. After fill tissue region with white and non-tissue with black in photoshop or other image processing software manually:
python Grid_filter.py

### 4. Analysis MISAR-seq data:
Rscript MISAR-seq.R

