#load packages
rm(list = ls())+gc()

pkgs <- c("Seurat","Signac","ArchR","SeuratData","SeuratWrappers","ggplot2","Matrix","stringr","grid","png","jpeg","export","patchwork","cowplot",
"harmony","ggalluvial","clusterProfiler","enrichplot","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene","DOSE","tidyverse","conflicted",
"ComplexHeatmap","Rsamtools","presto", "parallel","dplyr","BSgenome.Mmusculus.UCSC.mm10","TFBSTools","chromVARmotifs","motifmatchr")
suppressWarnings(suppressPackageStartupMessages(lapply(pkgs,library,character.only=TRUE)))

addArchRGenome("mm10")
addArchRThreads(threads = 20)


#Set species, sampleId and file_name
Species <- "mouse"
sampleId <- c("E11_0-S1","E13_5-S1","E15_5-S1","E18_5-S1")
file_name <- "Section1"

#Set color for each cluster
#For ATAC cluster color
Color_ATAC<-c("1"="#fcb462","2"="#9500ff","3"="#fb4f2b","4"="#4993b6","5"="#1e25ce","6"="#197a35","7"="#fd8000","8"="#d73027","9"="#878787","10"="#bb207c","11"="#ac6d03")
#For RNA and combined cluster color
Color_RNA<-c("1"="#fcb462","2"="#9500ff","3"="#fb4f2b","4"="#4993b6","5"="#1e25ce","6"="#197a35","7"="#fd8000","8"="#d73027","9"="#8e0054","10"="#bb207c","11"="#ac6d03",
"12"="#bc80bc","13"="#742c7e","14"="#05315e")
#show color
#scales::show_col(color_list)

#Set seed
set.seed(123)

#Set work directory
setwd("~/MISET-seq/Data")

#Read MISET-seq fragments file into ArchR
input_ATAC <- c()
for (i in sampleId){input_ATAC[i] <- paste0(i,'_filtered_fragments.tsv.bgz')}

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
  outputDirectory = file_name,
  #copyArrows = TRUE
)

#Add spatial location to proj
barcode <- read.csv("~/MISET-seq/barcode_file/MISET-seq_barcode.csv",header = T)
barcode_all<-rbind(barcode,barcode,barcode,barcode)
rownames(barcode_all) <- paste0(c(rep(sampleId[1],2500),rep(sampleId[2],2500),rep(sampleId[3],2500),rep(sampleId[4],2500)),"#",barcode_all$barcode)
barcode_final <- barcode_all[rownames(barcode_all)%in%rownames(proj),]
barcode_final<-barcode_final[,-1]

QC_index <- c("array_col","array_row")
for (i in QC_index) {
  proj<- addCellColData(proj, data = barcode_final[,i], cells = rownames(barcode_final), name = i, force = T)
}

#Add MISET-seq RNA into proj
#read RNA counts from h5 file
input_RNA <- c()
for (i in sampleId){input_RNA[i] <- paste0(i,'_raw_feature_bc_matrix.h5')}

seRNA <- import10xFeatureMatrix(
  input = input_RNA,
  names = sampleId
)
#Ignoring "Error in combining individual feature matrices! Returning as a list of individual feature matrices!"

#order of the features for each seRNA
conflict_prefer("rowRanges", "MatrixGenerics")
rowRanges(seRNA[[1]])<-rowRanges(seRNA[[2]])<-rowRanges(seRNA[[3]])<-rowRanges(seRNA[[4]])
seRNA_combined<-cbind(seRNA[[1]],seRNA[[2]],seRNA[[3]],seRNA[[4]])

#add RNA to ArchR object
proj <- addGeneExpressionMatrix(input = proj, seRNA = seRNA_combined)

#LSI-ATAC
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list(
    resolution = c(2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:40,
  force = TRUE,
  saveIterations = FALSE
)

#Remove ATAC batch with harmony
proj  <- addHarmony(
  ArchRProj = proj ,
  reducedDims = "IterativeLSI",
  name = "Harmony_ATAC",
  groupBy = "Sample"
)

# Find clusters for ATAC
proj <- addClusters(
  input = proj,
  reducedDims = "Harmony_ATAC", 
  method = "seurat",
  name = "Clusters_ATAC",
  resolution = 0.3,
  prefix = "",
  force = TRUE
)

#LSI-RNA
proj <- addIterativeLSI(
    ArchRProj = proj, 
    clusterParams = list(
      resolution = 2, 
      sampleCells = 10000,
      n.start = 10
    ),
    saveIterations = FALSE,
    useMatrix = "GeneExpressionMatrix", 
    depthCol = "Gex_nUMI",
    varFeatures = 2500,
    firstSelection = "variable",
    binarize = FALSE,
    name = "LSI_RNA"
)

#Remove RNA batch with harmony
proj  <- addHarmony(
  ArchRProj = proj ,
  reducedDims = "LSI_RNA",
  name = "Harmony_RNA",
  groupBy = "Sample",
  force = TRUE
)

# Find clusters for RNA
proj <- addClusters(
  input = proj,
  reducedDims = "Harmony_RNA", 
  method = "Seurat",
  name = "Clusters_RNA",
  resolution = 0.5,
  prefix = "",
  force = TRUE
)
}

#Combined Dims
proj <- addCombinedDims(proj, reducedDims = c("Harmony_ATAC", "Harmony_RNA"), name =  "Harmony_Combined")

#UMAPs
proj <- addUMAP(proj, reducedDims = "Harmony_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
proj <- addUMAP(proj, reducedDims = "Harmony_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)
proj <- addUMAP(proj, reducedDims = "Harmony_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)

# Find clusters of RNA+ATAC
proj <- addClusters(
  input = proj,
  reducedDims = "Harmony_Combined", 
  method = "Seurat",
  name = "Clusters_Combined",
  resolution = 0.5,
  prefix = "",
  force = TRUE
)


#save RNA+ ATAC ArchR project
saveRDS(proj,paste0(file_name,"_MISET-seq_ArchR.rds"))

#save sample in umap plot for each data
p1 <- plotEmbedding(proj, name = "Sample", embedding = "UMAP_RNA", rastr = FALSE, plotAs = "points", size = 0.5, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(proj, name = "Sample", embedding = "UMAP_ATAC", rastr = FALSE, plotAs = "points", size = 0.5, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(proj, name = "Sample", embedding = "UMAP_Combined",  rastr = FALSE, plotAs = "points", size = 0.5, labelAsFactors=F, labelMeans=F)
plotPDF(p1, p2, p3, name = "UMAP-RNA-ATAC-Combined-Sample", addDOC = FALSE)

p4 <- plotEmbedding(proj, name = c("Clusters_ATAC"), embedding = "UMAP_ATAC", rastr = FALSE, plotAs = "points", pal = color_list, size = 0.5, labelAsFactors=F, labelMeans=F)
p5 <- plotEmbedding(proj, name = c("Clusters_RNA"), embedding = "UMAP_RNA", rastr = FALSE, plotAs = "points", pal = color_list, size = 0.5, labelAsFactors=F, labelMeans=F)
p6 <- plotEmbedding(proj, name = c("Clusters_Combined"), embedding = "UMAP_Combined", rastr = FALSE, plotAs = "points", pal = color_list, size = 0.5, labelAsFactors=F, labelMeans=F)
plotPDF(p4, name = "UMAP-ATAC-Cluster", addDOC = FALSE) 
plotPDF(p5, name = "UMAP-RNA-Cluster", addDOC = FALSE) 
plotPDF(p6, name = "UMAP-Combined-Cluster", addDOC = FALSE)

#get each sample cellcoldata to plot Spatial cluster
idxSample <- list()
cellsSample <- list()
df <- list()
for(i in sampleId){
idxSample[[i]] <- BiocGenerics::which(proj$Sample %in% i)
cellsSample[[i]] <- proj$cellNames[idxSample[[i]]]
df[[i]] <- getCellColData(proj[cellsSample[[i]],])
}

pdf(file = paste0(file_name,"_Harmony_ATAC_Spatial_Cluster.pdf"), width=8.6, height=8.6)
for(i in sampleId){
imported_raster=OpenImageR::readImage(paste0("Data/", i, "-HE.jpg"))
g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)
p10 <- ggplot(as.data.frame(df[[i]]), aes(x = as.numeric(array_col), y = as.numeric(array_row), color=df[[i]][,"Clusters_ATAC"])) + scale_color_manual(values = Color_ATAC) +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    ggtitle(i) +
    geom_point(size = 2.5,shape=15) +
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
    scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
    coord_equal(xlim=c(0,51),ylim=c(51,1)) +
    theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
      axis.text=element_text(size=20),
      axis.title=element_text(size=20,face="bold"),
      legend.text=element_text(size=20),
      legend.title = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank())
  plot(p10)
}
dev.off()

pdf(file = paste0(file_name,"_Harmony_RNA_Spatial_Cluster.pdf"), width=8.6, height=8.6)
for(i in sampleId){
imported_raster=OpenImageR::readImage(paste0("Data/", i, "-HE.jpg"))
g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)
p7 <- ggplot(as.data.frame(df[[i]]), aes(x = as.numeric(array_col), y = as.numeric(array_row), color=df[[i]][,"Clusters_RNA"])) + scale_color_manual(values = Color_RNA) +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    ggtitle(i) +
    geom_point(size = 2.5, shape=15) +
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
    scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
    coord_equal(xlim=c(0,51),ylim=c(51,1)) +
    theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
      axis.text=element_text(size=20),
      axis.title=element_text(size=20,face="bold"),
      legend.text=element_text(size=20),
      legend.title = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank())
  plot(p7)
}
dev.off()

pdf(file = paste0(file_name,"_Harmony_Combined_Spatial_Cluster.pdf"), width=8.6, height=8.6)
for(i in sampleId){
imported_raster=OpenImageR::readImage(paste0("Data/", i, "-HE.jpg"))
g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)
p14 <- ggplot(as.data.frame(df[[i]]), aes(x = as.numeric(array_col), y = as.numeric(array_row), color=df[[i]][,"Clusters_Combined"])) + scale_color_manual(values = Color_RNA) +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    ggtitle(i) +
    geom_point(size = 2.5,shape=15) +
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
    scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
    coord_equal(xlim=c(0,51),ylim=c(51,1)) +
    theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
      axis.text=element_text(size=20),
      axis.title=element_text(size=20,face="bold"),
      legend.text=element_text(size=20),
      legend.title = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank())
  plot(p14)
}
dev.off()


#RNA and ATAC cluster frequency
pdf("ATAC_RNA_Cluster_frequency_heatmap.pdf", height=5.6, width=5.6)
A<-proj@cellColData$Clusters_RNA
B<-proj@cellColData$Clusters_ATAC
predictions <- table(A, B)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)
ggplot(predictions, aes(A, B, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
    low = "#ffffc8", high = "#7d0025") + xlab("Spatial-RNA-seq Clusters") + ylab("Spatial-ATAC-seq Clusters") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

#plot Sample Cluster overlap heatmap
pdf("Combined_Cluster_heatmap.pdf", height=5.6, width=5.6)
cM_Combined <- confusionMatrix(paste0(proj$Clusters_Combined), paste0(proj$Sample))
cM_Combined <- cM_Combined / Matrix::rowSums(cM_Combined)
cM_Combined <- as.matrix(cM_Combined)
p <- pheatmap::pheatmap(
  mat = cM_Combined[,order(colnames(cM_Combined))],
  cluster_rows = FALSE,
  cluster_cols = FALSE, 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)
p
dev.off()
