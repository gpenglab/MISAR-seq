## MISET-seq
### 1.Replace barcode in cellranger-arc
a. Enter cellranger-arc-2.0.1/lib/python/cellranger/barcodes/ and cellranger-arc-2.0.1/lib/python/atac/barcodes/;
b. Replace default barcode file "737K-arc-v1.txt.gz" with your new custom barcode file (in barcode_file, also named as "737K-arc-v1.txt.gz").

### 2.Generate fragments and raw_feature_bc matrix file for MISET-seq:
sh MISET-seq.sh

### 3.After fill tissue region with white and non-tissue with black in photoshop or other image processing software manually:
python Grid_filter.py

### 4.Analysis MISET-seq data in ArchR:
Rscripts MISET-seq-ArchR.R

### 5.Analysis MISET-seq data in Signac:
Rscripts MISET-seq-Signac.R
