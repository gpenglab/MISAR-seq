SelectGenes <- function(object,
                        atac.assay = "ATAC",
                        rna.assay = "RNA",
                        var.cutoff.gene = 0.9,
                        trajectory.name = "Trajectory",
                        distance.cutoff = 1000,
                        max.dist = 250000,
                        cor.cutoff = 0,
                        fdr.cutoff = 1e-04,
                        labelTop1 = 10,
                        labelTop2 = 10,
                        genome = "mm10") {
  trajRNA <- GetTrajectory(
    object,
    assay = rna.assay,
    slot = "data",
    trajectory.name = trajectory.name,
    smoothWindow = 7,
    log2Norm = TRUE
  )

  trajATAC <- GetTrajectory(
    object,
    assay = atac.assay,
    slot = "data",
    trajectory.name = trajectory.name,
    smoothWindow = 7,
    log2Norm = TRUE
  )

  # note here we only use the top 10% most variable genes
  groupMatRNA <- suppressMessages(
    TrajectoryHeatmap(
      trajRNA,
      varCutOff = var.cutoff.gene,
      pal = paletteContinuous(set = "horizonExtra"),
      limits = c(-2, 2),
      returnMatrix = TRUE
    )
  )

  groupMatATAC <- suppressMessages(
    TrajectoryHeatmap(
      trajATAC,
      varCutOff = NULL,
      maxFeatures = nrow(trajATAC),
      pal = paletteContinuous(set = "solarExtra"),
      limits = c(-2, 2),
      name = "Chromatin accessibility",
      returnMatrix = TRUE
    )
  )

  message("Linking cis-regulatory elements to genes...")
  df.p2g <- PeakToGene(peak.mat = groupMatATAC,
                       gene.mat = groupMatRNA,
                       genome = genome,
                       max.dist = max.dist)

  df.p2g <- df.p2g %>%
    subset(distance > distance.cutoff) %>%
    subset(Correlation > cor.cutoff & FDR < fdr.cutoff)

  trajATAC <- trajATAC[df.p2g$peak,]
  trajRNA <- trajRNA[df.p2g$gene,]

  return(res)

}

TrajectoryHeatmap <- function(trajectory,
                              varCutOff = 0.9,
                              maxFeatures = 25000,
                              scaleRows = TRUE,
                              rowOrder = NULL,
                              limits = c(-1.5, 1.5),
                              labelRows = FALSE,
                              pal = NULL,
                              labelMarkers = NULL,
                              labelTop = 50,
                              name = "Heatmap",
                              returnMatrix = FALSE) {
  mat <- assay(trajectory)

  #Rows with NA
  rSNA <- rowSums(is.na(mat))
  if (sum(rSNA > 0) > 0) {
    message("Removing rows with NA values...")
    mat <- mat[rSNA == 0, ]#Remove NA Rows
  }

  varQ <- ArchR:::.getQuantiles(matrixStats::rowVars(mat))

  orderedVar <- FALSE

  if (is.null(rowOrder)) {
    mat <- mat[order(varQ, decreasing = TRUE), ]
    orderedVar <- TRUE
    if (is.null(varCutOff) & is.null(maxFeatures)) {
      n <- nrow(mat)
    } else if (is.null(varCutOff)) {
      n <- maxFeatures
    } else if (is.null(maxFeatures)) {
      n <- (1 - varCutOff) * nrow(mat)
    } else{
      n <- min((1 - varCutOff) * nrow(mat), maxFeatures)
    }
    n <- min(n, nrow(mat))
    mat <- mat[head(seq_len(nrow(mat)), n),]
  }


  if (!is.null(labelTop) & labelTop > 0) {
    if (orderedVar) {
      idxLabel <- rownames(mat)[seq_len(labelTop)]
    } else{
      idxLabel <-
        rownames(mat)[order(varQ, decreasing = TRUE)][seq_len(labelTop)]
    }
  } else{
    idxLabel <- NULL
  }

  if (scaleRows) {
    mat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), `/`)
    mat[mat > max(limits)] <- max(limits)
    mat[mat < min(limits)] <- min(limits)
  }

  if (nrow(mat) == 0) {
    stop("No Features Remaining!")
  }

  if (is.null(pal)) {
    pal <- ArchR::paletteContinuous(set = "blueYellow", n = 100)
  }

  if (!is.null(rowOrder)) {
    idx <- rowOrder
  } else{
    idx <- order(unlist(apply(mat, 1, which.max)))
  }

  if(!is.null(idxLabel)){
    customRowLabel <- match(idxLabel, rownames(mat[idx,]))
  } else{
      customRowLabel <- NULL
  }

  ht <- ArchR:::.ArchRHeatmap(
    mat = mat[idx, ],
    scale = FALSE,
    limits = c(min(mat), max(mat)),
    color = pal,
    clusterCols = FALSE,
    clusterRows = FALSE,
    labelRows = labelRows,
    labelCols = FALSE,
    customRowLabel = customRowLabel,
    showColDendrogram = TRUE,
    name = name,
    draw = FALSE
  )

  if (returnMatrix) {
    return(mat[idx, ])
  } else{
    return(ht)
  }
}

PeakToGene <- function(peak.mat,
                       gene.mat,
                       genome = "hg19",
                       max.dist = 250000,
                       method = "corrleation") {
  if (!(genome %in% c("hg19", "hg38", "mm9", "mm10"))) {
    stop("Available genome are: hg19, hg38, mm9, and mm10!")
  }

  if (genome == "hg19") {
    gene_anno <- geneAnnoHg19
  } else if (genome == "hg38") {
    gene_anno <- geneAnnoHg38
  } else if (genome == "mm9") {
    gene_anno <- geneAnnoMm9
  } else if (genome == "mm10") {
    gene_anno <- geneAnnoMm10
  }

  ## create object for RNA data
  genes <- gene_anno$genes
  gene.use <-
    intersect(elementMetadata(genes)[, "symbol"], rownames(gene.mat))

  genes <- genes[elementMetadata(genes)[, "symbol"] %in% gene.use]
  gene.mat <- gene.mat[gene.use, ]

  gene_start <- ifelse(genes@strand == "+",
                       genes@ranges@start,
                       genes@ranges@start + genes@ranges@width - 1)

  genes <- GRanges(
    genes@seqnames,
    ranges = IRanges(gene_start,
                     width = 1),
    name = genes$symbol,
    gene_id = genes$gene_id,
    strand = genes@strand
  )

  seRNA <-
    SummarizedExperiment(assays = SimpleList(RNA = gene.mat),
                         rowRanges = genes)

  ## create object for ATAC data
  df_peak <-
    stringr::str_split_fixed(rownames(peak.mat), "-", 3)

  peakSet <- GRanges(df_peak[, 1],
                     IRanges(start = as.numeric(df_peak[, 2]),
                             end = as.numeric(df_peak[, 3])))

  seATAC <-
    SummarizedExperiment(assays = SimpleList(ATAC = peak.mat),
                         rowRanges = peakSet)

  ## find putative peak-to-gene
  o <-
    data.frame(findOverlaps(
      resize(seRNA, 2 * max.dist + 1, "center"),
      resize(rowRanges(seATAC), 1, "center"),
      ignore.strand = TRUE
    ))
  o$distance <- IRanges::distance(rowRanges(seRNA)[o[, 1]],
                                  rowRanges(seATAC)[o[, 2]])
  colnames(o) <- c("gene_idx", "peak_idx", "distance")

  df <- rowRanges(seATAC)[o$peak_idx, ]

  o$gene <- rowData(seRNA)[o$gene_idx, ]$name
  o$peak <- paste0(
    df@seqnames,
    "-",
    as.data.frame(df@ranges)$start,
    "-",
    as.data.frame(df@ranges)$end
  )


  ## compute correlation
  o$Correlation <- rowCorCpp(as.integer(o$peak_idx),
                             as.integer(o$gene_idx),
                             assay(seATAC),
                             assay(seRNA))

  ## compute p-value
  o$TStat <-
    (o$Correlation / sqrt((
      pmax(1 - o$Correlation ^ 2, 0.00000000000000001, na.rm = TRUE)
    ) / (ncol(seATAC) - 2))) #T-statistic P-value

  o$Pval <- 2 * pt(-abs(o$TStat), ncol(seATAC) - 2)
  o$FDR <- p.adjust(o$Pval, method = "fdr")
  o <- o[!is.na(o$FDR), ]

  return(o)

}

GetTrajectory <- function(object = NULL,
                          trajectory.name = "Trajectory",
                          assay = NULL,
                          slot = "counts",
                          groupEvery = 1,
                          log2Norm = TRUE,
                          scaleTo = 10000,
                          smoothWindow = 11) {
  if (is.null(assay) | !assay %in% Assays(object)) {
    stop("Please provide an available assay!")
  }

  if (!(trajectory.name %in% colnames(object@meta.data))) {
    stop(glue::glue("Cannot find trajecotry {trajectory.name}!"))
  }

  trajectory <- object@meta.data[trajectory.name]

  trajectory <- trajectory[!is.na(trajectory[, 1]), , drop = FALSE]
  breaks <- seq(0, 100, groupEvery)
  if (!all(is.numeric(trajectory[, 1]))) {
    stop("Trajectory must be a numeric. Did you add the trajectory with addTrajectory?")
  }
  if (!all(trajectory[, 1] >= 0 & trajectory[, 1] <= 100)) {
    stop(
      "Trajectory values must be between 0 and 100. Did you add the trajectory with addTrajectory?"
    )
  }

  groupList <- lapply(seq_along(breaks), function(x) {
    if (x == 1) {
      NULL
    } else{
      rownames(trajectory)[which(trajectory[, 1] > breaks[x - 1] &
                                   trajectory[, 1] <= breaks[x])]
    }
  })[-1]

  names(groupList) <-
    paste0("T.", breaks[-length(breaks)], "_", breaks[-1])

  message("Creating Trajectory Group Matrix..")
  data.use <- GetAssayData(object, assay = assay, slot = slot)

  groupMat <- lapply(1:length(groupList), function(x) {
    cell_names <- groupList[[x]]
    mat <- Matrix::rowMeans(data.use[, cell_names])

  }) %>% Reduce(cbind, .)

  colnames(groupMat) <- names(groupList)

  #Scale
  if (!is.null(scaleTo)) {
    if (any(groupMat < 0)) {
      message(
        "Some values are below 0, this could be the Motif activity matrix in which scaleTo should be set = NULL.\nContinuing without depth normalization!"
      )
    } else{
      groupMat <- t(t(groupMat) / colSums(groupMat)) * scaleTo
    }
  }

  if (log2Norm) {
    if (any(groupMat < 0)) {
      message(
        "Some values are below 0, this could be a Motif activity matrix in which log2Norm should be set = FALSE.\nContinuing without log2 normalization!"
      )
    } else{
      groupMat <- log2(groupMat + 1)
    }
  }

  if (!is.null(smoothWindow)) {
    message("Smoothing...")
    smoothGroupMat <-
      as.matrix(t(apply(groupMat, 1, function(x)
        centerRollMean(x, k = smoothWindow))))
    colnames(smoothGroupMat) <- paste0(colnames(groupMat))
    colnames(groupMat) <- paste0(colnames(groupMat))

    #Create SE
    seTrajectory <- SummarizedExperiment(assays = SimpleList(
      smoothMat = as.matrix(smoothGroupMat),
      mat = as.matrix(groupMat)
    ))

  } else{
    colnames(groupMat) <- paste0(colnames(groupMat))

    #Create SE
    seTrajectory <-
      SummarizedExperiment(assays = SimpleList(mat = as.matrix(groupMat)))

  }

  return(seTrajectory)
}

getScores <- function(ATAC.se,
                      DERTab,
                      geneList,
                      nCores=4){

  DERTab <- DERTab[DERTab$Gene %in% geneList,] # Filter only to these genes
  DERGenes <- sort(as.character(unique(DERTab$Gene)))

  ATAC.mat <- assay(ATAC.se)

  time_elapsed <- Sys.time()

  DERMatL <- pbmcapply::pbmclapply(X=DERGenes,
                                    FUN=function(x) {

                                      DERPeaks <- unique(DERTab$Peak[DERTab$Gene %in% x])

                                      if(length(DERPeaks) > 1) {
                                        DERCounts <- Matrix::colSums(ATAC.mat[DERPeaks,])
                                      } else if(length(DERPeaks==1)) {
                                        DERCounts <- ATAC.mat[DERPeaks,]
                                      }
                                    },mc.cores = nCores)

  DERMat <- Matrix::Matrix(do.call('rbind',DERMatL),sparse=TRUE)

  rownames(DERMat) <- DERGenes

  return(DERMat)

}


### prepare data
cortex <- readRDS("data/seurat_rna_atac_all.rds")

DefaultAssay(cortex) <- "peaks"
sparse_mat <- cortex@assays$peaks@counts
rowidx <- which(rowSums(sparse_mat > 0) > 5)
sparse_mat <- sparse_mat[rowidx, ]
toexportpeaks <- cortex@assays$peaks@ranges[rowidx]
names(toexportpeaks) = rownames(sparse_mat)
fragment_counts <- SummarizedExperiment(assays = list(counts = sparse_mat), rowRanges = toexportpeaks)
assayNames(fragment_counts) <- "counts"
ATAC.se <- fragment_counts

DefaultAssay(cortex) <- "RNA"
cortex <- NormalizeData(cortex)
rna_matrix <- cortex@assays$RNA@data
rowidx <- which(rowSums(rna_matrix > 0) > 5)
rna_matrix <- rna_matrix[rowidx, ]
RNAmat <- rna_matrix


### get DER score
res <- SelectGenes(object = cortex,
                        atac.assay = "peaks",
                        rna.assay = "RNA",
                        var.cutoff.gene = 0,
                        trajectory.name = "Pseudotime_Combined")
p2g <- res
p2g.cis = data.frame(Peak = match(p2g$peak, rownames(sparse_mat)), PeakRanges = p2g$peak, Gene = p2g$gene, rObs = p2g$Correlation, pvalZ = p2g$Pval)
p2g.cis = na.omit(p2g.cis)

cisCorr <- p2g.cis
cisCorr = p2g.cis[grep("*Rik|^Gm|^mt-|^Rps|^Rpl", p2g.cis$Gene, invert = T),]
cisCorr <- cisCorr %>% group_by(Peak) %>% filter(rObs==max(rObs))
cisCorr.filt <- cisCorr %>% filter(pvalZ <= 0.05)

DERMat <- getScores(
    ATAC.se = ATAC.se,
    DERTab = cisCorr.filt,
    geneList = DERGenes,
    nCores = 10
)



