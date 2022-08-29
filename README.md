## MISET-seq
# Replace barcode in cellranger-arc
1. Enter cellranger-arc-2.0.1/lib/python/cellranger/barcodes/ and cellranger-arc-2.0.1/lib/python/atac/barcodes/
2. Replace default barcode file "737K-arc-v1.txt.gz" with your new custom barcode file (in barcode_file, also named as "737K-arc-v1.txt.gz").

# Generate fragments and raw_feature_bc matrix file for MISET-seq:
sh MISET-seq.sh

# After fill tissue region with white and non-tissue with black in photoshop or other image processing software manually:
python Grid_filter.py

# Analysis MISET-seq data in ArchR:
Rscripts MISET-seq-ArchR.R

# Analysis MISET-seq data in Signac:
Rscripts MISET-seq-Signac.R
