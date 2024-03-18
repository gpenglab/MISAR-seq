## MISAR-seq
### The first high-throughput spatial co-omics method, named as microfluidic indexing based spatial ATAC and RNA sequencing (MISAR-seq).
All RNA and ATAC fastq data can be download from NCBI: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP491963&o=acc_s%3Aa

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

PS1: New version in ZENODO: https://doi.org/10.5281/zenodo.7480069

PS2: All section1 section2 and E15_5-S1 data can download from Google Drive: https://drive.google.com/drive/folders/1xIUkhB-W2O3wHca0yfMaKAVXdVTBuf5N?usp=drive_link
