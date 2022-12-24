pkgs <- c("Seurat","Signac","ArchR","tidyverse","Matrix","patchwork",
          "ggalluvial","clusterProfiler","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene",
          "parallel","BSgenome.Mmusculus.UCSC.mm10","TFBSTools","chromVARmotifs","motifmatchr",
          "magrittr","viridis","igraph", "clustree", "reticulate")
suppressWarnings(suppressPackageStartupMessages(lapply(pkgs,library,character.only=TRUE)))

# set work directory
setwd("~/MISAR-seq")

# source function needed
source_python("Python/utils.py")
source("R/function.R")

sampleId <- "E15_5-S1"
color_list<-c("1"="#644498","2"="#E94F2F","3"="#488CAD","4"="#D6ACCF","5"="#207639",
              "6"="#EF7D18","7"="#7184C1","8"="#70B1D7","9"="#DCBE70","10"="#A66C22",
              "11"="#1A7F85","12"="#ED7C7A","13"="#A8CD92","14"="#A91D30","15"="#F1CC32",
              "16"="#E6E754","17"="#063D20","18"="#8dd3c8","19"="#b31631","20"="#fbd326"
              )

set.seed(123)
addArchRGenome("mm10")
addArchRThreads(threads = 30)

# read counts
counts <- Read10X_h5("Data/E15_5-S1_raw_feature_bc_matrix.h5")

# create a Seurat object
dat <- CreateSeuratObject(
  project = sampleId,
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# calculate percent of mitochondrial for RNA
dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^mt-")

# read position file for filter
location <- read.table(paste0("Data/position_", sampleId,".txt"), sep =",", header = FALSE, dec =".", stringsAsFactors = F)
x <- as.character(location[1,])[-1]
x <- as.data.frame(x)
colnames(x) <- "array"
y <- read.csv("Barcode/MISAR-seq_barcode_filter.csv", header = T) 
z <- merge(x,y,by="array")

# filter non-tissue barcodes
dat <- subset(dat, cells=z$barcode)
dat <- SCTransform(dat, vars.to.regress = "percent.mt", do.scale = T, verbose = FALSE) %>% RunPCA(verbose = FALSE)

# read fragments matrix and filter non-tissue barcodes
frag <- read.table("Data/atac_fragments.tsv.gz", sep = "\t")
colnames(frag) <- c("chr","start","end","barcode","fragments")
frag <- frag[frag$barcode %in% z$barcode,]

write.table(frag, file = paste0("Data/",sampleId,'_filtered_fragments.tsv'), sep = '\t',col.names=FALSE, row.names = FALSE,quote = FALSE)
rm(frag)
# compress tabix file with bgzip for index
bgzip(paste0("Data/",sampleId,'_filtered_fragments.tsv'))
file.remove(paste0("Data/",sampleId,"_filtered_fragments.tsv"))

# create ArchR object
input_ATAC <- paste0("Data/",sampleId,'_filtered_fragments.tsv.bgz')

ArrowFiles <- createArrowFiles(
  inputFiles = input_ATAC,
  sampleNames = sampleId,
  minTSS = 0,
  minFrags = 0,
  maxFrags = 1e+07,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  offsetPlus = 0,
  offsetMinus = 0,
  #force = TRUE,
  TileMatParams = list(tileSize = 5000)
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = sampleId,
  #copyArrows = TRUE
)

#read RNA counts from h5 file
input_RNA <- paste0("Data/",sampleId,'_raw_feature_bc_matrix.h5')

seRNA <- import10xFeatureMatrix(
  input = input_RNA,
  names = sampleId,
  featureType = "Gene Expression" # or "Peaks" for ATAC
)

# add RNA data to proj
proj <- addGeneExpressionMatrix(input = proj, seRNA = seRNA)

# read the spatial barcode
barcode <- read.csv("Barcode/MISAR-seq_barcode.csv",header = T,row.names = 1)

# filter non-tissue barcode
barcode <- barcode[z$barcode,]

# read and filter per_barcode_metrics.csv file which from cellranger
metrics <- read.csv(paste0("Data/",sampleId,"_per_barcode_metrics.csv"), header = T, row.names = 1)
metrics <- merge(metrics,barcode,by="row.names")
rownames(metrics) <- metrics$Row.names
rownames(metrics) <- paste0(sampleId,"#",rownames(metrics))

metrics <- transform(metrics, percent_mito_atac=(metrics$atac_mitochondrial_reads/metrics$atac_raw_reads)*100, 
                     FRiP=(metrics$atac_peak_region_fragments/metrics$atac_fragments)*100, 
                     percent_TSS_fragments=(metrics$atac_TSS_fragments/metrics$atac_fragments)*100, 
                     total_fragments=metrics$atac_fragments, Sample=sampleId)

QC_index <- c("array_col","array_row","percent_mito_atac","FRiP","percent_TSS_fragments","total_fragments")

# add QC_index into proj
for (i in QC_index) {
  proj <- addCellColData(ArchRProj = proj, data = metrics[,i], cells = rownames(metrics),name = i, force = T)
}

# plot RNA and ATAC spatial feature
pdf(file = paste0(sampleId,"_Spatial_QC.pdf"), width=8.6, height=8.6)

df1 <- as.data.frame(proj@cellColData)
df1 <- transform(df1, log10_unique_fragments=log10(df1$nFrags))

imported_raster=OpenImageR::readImage(paste0("Data/",sampleId,"-HE.jpg"))
g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)

for (i in c("Gex_nUMI","Gex_nGenes","Gex_MitoRatio","Gex_RiboRatio","TSSEnrichment","nFrags","log10_unique_fragments","total_fragments","FRiP","percent_mito_atac","percent_TSS_fragments")){
  p1 <- ggplot(df1, aes(x = as.numeric(array_col), y = as.numeric(array_row), color=df1[,i])) + scale_color_gradientn(colours = ArchRPalettes$solarExtra) +
    ggtitle(ylab) +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    ggtitle(i) +
    geom_point(size = 2.5, shape=15) +
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(name=NULL, limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
    scale_y_reverse(name=NULL, limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
    coord_equal(xlim=c(0,51),ylim=c(51,1)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          axis.text=element_text(size=0),
          axis.ticks = element_blank(),
          legend.text=element_text(size=20),
          legend.title = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
  plot(p1)
}
dev.off()

# using spatial position to smooth
# for ATAC
proj <- IterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix", 
  name = "LSI_ATAC", 
  iterations = 2, 
  clusterParams = list(
    resolution = c(2),
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force = TRUE,
  saveIterations = FALSE,
  verbose = F
)

# for RNA
proj <- IterativeLSI(
  ArchRProj = proj, 
  iterations = 3, 
  clusterParams = list(
    resolution = 2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  varFeatures = 3000,
  dimsToUse = 1:20,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA",
  outlierQuantiles = c(0.02, 0.98),
  filterQuantile = 0.995,
  force = T, 
  verbose = F
)

# repalce ArchR RNA reducedDims with Seurat
matSVD <- dat@reductions$pca@cell.embeddings[,1:30]
colnames(matSVD) <- gsub("PC_", "LSI", colnames(matSVD))
rownames(matSVD) <- paste0(sampleId,"#",rownames(matSVD))
matSVD <- matSVD[rownames(proj),]
proj@reducedDims$LSI_RNA$matSVD <- matSVD

# combine RNA and ATAC
proj <- addCombinedDims(proj, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")

# add UMAP for RNA, ATAC and Combined
for (i in c("ATAC","RNA","Combined")) {
proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = paste0("LSI_", i),
  name = paste0("UMAP_",i),
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE)
}

# perform clustering for RNA, ATAC and Combined
for (i in c(0.3, 0.5, 0.8)){
  proj <-  IterativeLSI_Clustering(
    input = proj,
    reducedDims = "LSI_RNA", 
    method = "Seurat",
    name = paste0("Clusters_RNA_",i),
    resolution = i,
    prefix = "",
    force = TRUE,
    biasCol = "Gex_nUMI",
    knnAssign = 30,
    nOutlier = 20,
    filterBias = T,
    verbose = F)

  proj <-  IterativeLSI_Clustering(
    input = proj,
    reducedDims = "LSI_ATAC", 
    method = "Seurat",
    name = paste0("Clusters_ATAC_",i),
    resolution = i,
    prefix = "",
    force = TRUE,
    knnAssign = 30,
    nOutlier = 20,
    filterBias = T,
    verbose = F)

  proj <-  IterativeLSI_Clustering(
    input = proj,
    reducedDims = "LSI_Combined", 
    method = "Seurat",
    name = paste0("Clusters_Combined_",i),
    resolution = i,
    prefix = "",
    knnAssign = 30,
    nOutlier = 20,
    filterBias = T,
    force = TRUE,
    verbose = F
  )
}

# plot clusters in UMAP for RNA, ATAC and Combined
df2 <- proj@cellColData

p2 <- plotEmbedding(proj, name = c("Clusters_ATAC_0.5"), embedding = "UMAP_ATAC", 
                    rastr = FALSE, plotAs = "points", 
                    pal = color_list[sort(as.numeric(unique(df2[,"Clusters_ATAC_0.5"])))],
                    size = 0.3, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(proj, name = c("Clusters_RNA_0.5"), embedding = "UMAP_RNA",
                    rastr = FALSE, plotAs = "points", 
                    pal = color_list[sort(as.numeric(unique(df2[,"Clusters_RNA_0.5"])))], 
                    size = 0.3, labelAsFactors=F, labelMeans=F)
p4 <- plotEmbedding(proj, name = c("Clusters_Combined_0.5"), embedding = "UMAP_Combined",
                    rastr = FALSE, plotAs = "points", 
                    pal = color_list[sort(as.numeric(unique(df2[,"Clusters_Combined_0.5"])))], 
                    size = 0.3, labelAsFactors=F, labelMeans=F)

plotPDF(p2, name = "UMAP-ATAC-Cluster_0.5", addDOC = FALSE) 
plotPDF(p3, name = "UMAP-RNA-Cluster_0.5", addDOC = FALSE) 
plotPDF(p4, name = "UMAP-Combined-Cluster_0.5", addDOC = FALSE)

# plot spatial clusters for RNA, ATAC and Combined
pdf(file = paste0(sampleId,"_Spatial_Cluster_res_0.5.pdf"), width=8.6, height=8.6)

for(i in c("Clusters_ATAC_0.5","Clusters_RNA_0.5","Clusters_Combined_0.5")){
  imported_raster=OpenImageR::readImage(paste0("Data/", sampleId, "-HE.jpg"))
  g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)
  p5 <- ggplot(as.data.frame(df2), aes(x = as.numeric(array_col), y = as.numeric(array_row),
                                           color=df2[,i])) + 
    scale_color_manual(values = color_list[sort(as.numeric(unique(df2[,i])))],
                       labels=sort(as.numeric(unique(df2[,i])))) +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    ggtitle(i) +
    geom_point(size = 2.5, shape=15) +
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(name=NULL, limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
    scale_y_reverse(name=NULL, limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
    coord_equal(xlim=c(0,51),ylim=c(51,1)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          axis.text=element_text(size=0),
          axis.ticks = element_blank(),
          legend.text=element_text(size=20),
          legend.title = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
  plot(p5)
}
dev.off()
