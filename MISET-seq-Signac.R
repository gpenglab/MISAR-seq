
pkgs <- c("Signac","Seurat","Rsamtools","GenomeInfoDb","EnsDb.Hsapiens.v75","EnsDb.Mmusculus.v79","BSgenome.Hsapiens.UCSC.hg38","BSgenome.Mmusculus.UCSC.mm10","future","ggplot2","dplyr","reshape2","stringr")
suppressWarnings(suppressPackageStartupMessages(lapply(pkgs,library,character.only=TRUE)))

#Set sampleId
sampleId <- "E15_5-S1"

#Set seed
set.seed(123)

# load the RNA and ATAC data
counts <- Read10X_h5("raw_feature_bc_matrix.h5")
fragpath <- "atac_fragments.tsv.gz"

# Get gene annotations for mm10
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"

# create a Seurat object containing the RNA adata
df <- CreateSeuratObject(
  project = sampleId,
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
df[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks, sep = c(":", "-"), fragments = fragpath, annotation = annotation)

DefaultAssay(df) <- "RNA"

# Remove all mitochondrial gene
# df <- df[!grepl("^mt-", rownames(df)), ]

#Caculate the percent of mitochondrial and ribosome protein gene
df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^mt-")
df[["percent.rp"]] <- PercentageFeatureSet(df, pattern = "^Rp[sl]")

DefaultAssay(df) <- "ATAC"

df <- NucleosomeSignal(df)

#setting fast=TRUE will only compute the TSS enrichment score without storing the entire cell by position matrix of Tn5 insertion frequency for each grid
df <- TSSEnrichment(df, fast = FALSE)

# The set of peaks identified using Cellranger often merges distinct peaks that are close together. This can create a problem for certain analyses, particularly motif enrichment analysis and peak-to-gene linkage.
# call peaks using MACS2
peaks <- CallPeaks(df, macs2.path = "/MISET-seq/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions, mouse:blacklist_mm10, human:blacklist_hg38_unified.
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")

peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(df),
  features = peaks,
  cells = colnames(df)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
df[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

# read spatial barcode file.
barcode <- read.csv("/MISET-seq/barcode_file/MISET-seq_barcode.csv", header = T, row.names = 1)

# read and filter per_barcode_metrics.csv file which from cellranger
metrics <- read.csv("per_barcode_metrics.csv", header = T, row.names = 1)

metrics <- merge(metrics,barcode,by="row.names")

rownames(metrics) <- metrics$Row.names

metrics <- transform(metrics, atac_percent_mito = (metrics$atac_mitochondrial_reads/metrics$atac_raw_reads)*100, FRiP = (metrics$atac_peak_region_fragments/metrics$atac_fragments)*100, percent_TSS_fragments = (metrics$atac_TSS_fragments/metrics$atac_fragments)*100, log10_unique_fragments = log10(metrics$atac_fragments), Sample =sampleId)

metrics <- metrics[,c("FRiP", "log10_unique_fragments", "atac_percent_mito", "percent_TSS_fragments", "Sample", "array_col", "array_row" )]

df <- AddMetaData(object = df, metadata = metrics)

# filter non-tissue grids use filtered spatial position file
location <- read.table(paste0("position_", sampleId,".txt"), sep =",", header = FALSE, dec =".", stringsAsFactors = F)
x <- as.character(location[1,])[-1]
x <- as.data.frame(x)
colnames(x) <- "array"
y <- read.csv("/MISET-seq/barcode_file/MISET-seq_barcode_filter.csv", header = T) 
z <- merge(x,y,by="array")

df <- subset(df, cells=z$barcode)

#Plot QC data
pdf(paste0(sampleId, "_QC.pdf"), width=30, height=6)

VlnPlot(
  object = df,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rp", "nCount_ATAC","nCount_peaks", "TSS.enrichment", "nucleosome_signal", "FRiP", "log10_unique_fragments", "atac_percent_mito", "percent_TSS_fragments"),
  ncol = 12,
  pt.size = 0
)

dev.off()

#Plot QC data in spatial

#import tissue figure for spatial plot
imported_raster=OpenImageR::readImage(paste0(sampleId,"-HE.jpg"))
g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)

pdf(file = paste0(sampleId,"_Spatial_Feature.pdf"), width=8.6, height=8.6)

for (i in c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rp", "nCount_ATAC", "nCount_peaks", "TSS.enrichment", "nucleosome_signal", "FRiP", "log10_unique_fragments", "atac_percent_mito", "percent_TSS_fragments")){
  
  xlab=mean(df@meta.data[,i])
  
  p1 <- ggplot(df@meta.data,aes(x=orig.ident,y=df@meta.data[,i],fill=orig.ident))+geom_violin(alpha=0.8,width=1) + guides(fill=F)+xlab(paste0("mean=",round(xlab,2)))+ylab(i)+theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"), axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"), panel.background = element_blank()) + geom_boxplot(width=.2,col="black",fill="white") + NoLegend()
  
  p2 <- ggplot(df@meta.data, aes(x = as.numeric(array_col), y = as.numeric(array_row), color=df@meta.data[,i])) + scale_color_gradientn(colours = c("white", "red")) +
    ggtitle(i) +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    geom_point(size = 0.7, shape=15)+
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
    scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
    coord_equal(xlim=c(0,51),ylim=c(51,1)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20,face="bold"),
          legend.text=element_text(size=20),
          legend.title = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
  plot(p1 | p2) 
}

dev.off()

#perform RNA reduction and cluster analysis
DefaultAssay(df) <- "RNA"
df <- SCTransform(df, vars.to.regress = "percent.mt", verbose = T) %>% RunPCA(verbose = FALSE) %>% FindNeighbors(dims = 1:15) %>% FindClusters(resolution = 0.5) %>% RunUMAP(dims = 1:15, reduction.name = "umap.rna")

#perform ATAC reduction and cluster analysis
#The set of peaks identified using Cellranger often merges distinct peaks that are close together. This can create a problem for certain analyses, particularly motif enrichment analysis and peak-to-gene linkage. So we use macs2 to identify a more accurate set of peaks.
DefaultAssay(df) <- "peaks"
df <- FindTopFeatures(df, min.cutoff = 5)
df <- RunTFIDF(df)
df <- RunSVD(df)

#The first LSI component often captures sequencing depth (technical variation) rather than biological variation. 
#DepthCor(df, n=30)

df <- FindNeighbors(df, reduction = 'lsi', dims = 2:30) %>% FindClusters(resolution = 0.5, verbose = FALSE, algorithm = 3) %>% RunUMAP(reduction = 'lsi', dims = 2:30, reduction.name = "umap.atac")

#plot Cluster in each resolution
pdf(file = paste0(sampleId,"_Cluster_each_raw.pdf"), width=9.0, height=4.0)

p3 <- DimPlot(df, reduction ="umap.rna", group.by = "SCT_snn_res.0.5", label = T, pt.size = 1.5,label.size=4) + ggtitle(paste0("scRNA_seq_",i)) + theme(plot.title = element_text(hjust = 0.5,size = 15))

p4 <- DimPlot(df, reduction ="umap.atac", group.by ="peaks_snn_res.0.5", label = T, pt.size = 1.5, label.size=4) + ggtitle(paste0("scATAC_seq_",i)) + theme(plot.title = element_text(hjust = 0.5,size = 15))

plot(p3 | p4)

dev.off()


pdf(file = paste0(sampleId,"_Spatial_Clusters.pdf"), width=8.6, height=8.6)

for (i in c("SCT_snn_res.0.5","peaks_snn_res.0.5")){
  
  p5 <- ggplot(df@meta.data, aes(x = as.numeric(array_col), y = as.numeric(array_row), color=df@meta.data[,i])) +
    ggtitle(i) +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    geom_point(size = 2.5,shape=15)+
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
    scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
    coord_equal(xlim=c(0,51),ylim=c(51,1)) +
    theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20,face="bold"),
          legend.text=element_text(size=20),
          legend.title = element_blank(),
          #legend.title = element_text(colour="black", size=15, face="bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
  plot(p5)
}

dev.off()


#save seurat object
saveRDS(df, file = paste0(sampleId,"_Seurat.rds"))

