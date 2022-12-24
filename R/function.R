IterativeLSI <- function(
  ArchRProj = NULL, 
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list(
      resolution = c(2), 
      sampleCells = 10000,
      # maxClusters = 6,
      n.start = 10
  ),
  firstSelection = "top",
  depthCol = "nFrags",
  varFeatures = 25000,
  dimsToUse = 1:30,
  LSIMethod = 2,
  scaleDims = TRUE,
  corCutOff = 0.75,
  binarize = TRUE,
  outlierQuantiles = c(0.02, 0.98),
  filterBias = TRUE,
  sampleCellsPre = 10000,
  projectCellsPre = FALSE,
  sampleCellsFinal = NULL,
  selectionMethod = "var",
  scaleTo = 10000,
  totalFeatures = 500000,
  filterQuantile = 0.995,
  excludeChr = c(),
  saveIterations = TRUE,
  UMAPParams = list(
    n_neighbors = 40, 
    min_dist = 0.4, 
    metric = "cosine", 
    verbose = FALSE, 
    fast_sgd = TRUE
  ),
  nPlot = 10000,
  outDir = getOutputDirectory(ArchRProj),
  threads = getArchRThreads(),
  seed = 1,
  verbose = TRUE,
  force = FALSE,
  logFile = NULL,
  histology = FALSE,
  refine = "V2",
  refine_iter = 1
  ){


  if(varFeatures < 1000){
    stop("Please provide more than 1000 varFeatures!")
  }


  ArchR:::.requirePackage("Matrix", source = "cran")
  tstart <- Sys.time()

  if(!is.null(ArchRProj@reducedDims[[name]])){
    if(!force){
      stop("Error name in reducedDims Already Exists! Set force = TRUE or pick a different name!")
    }
  }

  #Set Seed
  set.seed(seed)
  outDir <- file.path(outDir, name)
  dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

  #All the Cell Names
  cellNames <- rownames(getCellColData(ArchRProj))
  if(!is.null(sampleCellsPre)){
    if(length(cellNames) < sampleCellsPre){
      sampleCellsPre <- NULL
    }
  }
  if(!is.null(sampleCellsFinal)){
    if(length(cellNames) < sampleCellsFinal){
      sampleCellsFinal <- NULL
    }
  }

  #MatrixFiles
  ArrowFiles <- getSampleColData(ArchRProj)[,"ArrowFiles"]

  #Check if Matrix is supported and check type
  if(tolower(useMatrix) == "tilematrix"){
    useMatrix <- "TileMatrix"
    tileSizes <- lapply(ArrowFiles, function(x){
      h5read(x, "TileMatrix/Info/Params/")$tileSize[1]
    }) %>% unlist
    if(length(unique(tileSizes)) != 1){
      stop("Error not all TileMatrices are the same tileSize!")
    }
    tileSize <- unique(tileSizes)
  }else if(tolower(useMatrix) == "peakmatrix"){
    useMatrix <- "PeakMatrix"
    tileSize <- NA
  }else{
    tileSize <- NA
  }

  units <- unique(unlist(lapply(ArrowFiles, function(x) h5read(x, paste0(useMatrix, "/Info/Units")))))
  if(length(units) != 1){
    stop("Units of matrices are not identical!")
  }
  if(grepl("log",units,ignore.case=TRUE)){
    stop("Cannot use log transformed values for iterativeLSI!")
  }

  tstart <- Sys.time()  
  ############################################################################################################################
  # Organize Information for LSI
  ############################################################################################################################
  chrToRun <- ArchR:::.availableSeqnames(ArrowFiles, subGroup = useMatrix)

  if(tolower(firstSelection) == "top"){
    
    if(!binarize){
      matClass <- h5read(ArrowFiles[1], paste0(useMatrix,"/Info/Class"))
      if(matClass != "Sparse.Binary.Matrix"){
        stop("Input matrix is not binarized and binarize != TRUE. Please use binarized data if using top selection for first iteration! Set binarize = TRUE!")
      }
    }

    #Compute Row Sums Across All Samples
    message(sprintf("%s : %s", Sys.time(), "Computing Total Across All Features"))
    if(useMatrix == "TileMatrix"){
      totalAcc <- ArchR:::.getRowSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun, addInfo = FALSE)
      totalAcc$start <- (totalAcc$idx - 1) * tileSize
    }else{
      totalAcc <- ArchR:::.getRowSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun, addInfo = TRUE)
    }

    #Filter Chromosomes
    if(length(excludeChr) > 0){
      totalAcc <- totalAcc[BiocGenerics::which(totalAcc$seqnames %bcni% excludeChr), , drop = FALSE]
    }

    #Identify the top features to be used here
    message(sprintf("%s : %s", Sys.time(), "Computing Top Features"))
    nFeature <- varFeatures[1]
    rmTop <- floor((1-filterQuantile) * totalFeatures)
    topIdx <- head(order(totalAcc$rowSums, decreasing=TRUE), nFeature + rmTop)[-seq_len(rmTop)]
    topFeatures <- totalAcc[sort(topIdx),]

    gc()

  }else if(tolower(firstSelection) %in% c("var", "variable")){

    if(binarize){
      stop("Please do not binarize data if using variable selection for first iteration! Set binarize = FALSE!")
    }

    if(units %in% "BinarizedCounts"){
      stop("Cannot do variable selection with BinarizedCounts. Set firstSelection = Top!")
    }

    #Compute Row Sums Across All Samples
    message(sprintf("%s : %s", Sys.time(), "Computing Variability Across All Features"))
    if(useMatrix == "TileMatrix"){
      totalAcc <- ArchR:::.getRowVars(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun, useLog2 = TRUE)
      totalAcc$start <- (totalAcc$idx - 1) * tileSize
    }else{
      totalAcc <- ArchR:::.getRowVars(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun, useLog2 = TRUE)
    }

    #Filter Chromosomes
    if(length(excludeChr) > 0){
      totalAcc <- totalAcc[BiocGenerics::which(totalAcc$seqnames %bcni% excludeChr), , drop = FALSE]
    }

    #Identify the top features to be used here
    message(sprintf("%s : %s", Sys.time(), "Computing Variable Features"))
    nFeature <- varFeatures[1]
    if(nFeature > 0.5 * nrow(totalAcc)){
      stop("nFeature for variable selection must be less than 1/2 the total features!")
    }
    topIdx <- head(order(totalAcc$combinedVars, decreasing=TRUE), nFeature)
    topFeatures <- totalAcc[sort(topIdx),]

    gc()

  }else{

    stop("firstSelect method must be Top or Var/Variable!")

  }

  cellDepth <- tryCatch({
      df <- getCellColData(ArchRProj = ArchRProj, select = depthCol)
      v <- df[,1]
      names(v) <- rownames(df)
      v
    }, error = function(e){
      tryCatch({
        ArchR:::.getColSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun)[ArchRProj$cellNames]
      }, error = function(y){
        stop("Could not determine depth from depthCol or colSums!")
      })
    }
  )
  cellDepth <- log10(cellDepth + 1)

  ############################################################################################################################
  # LSI Iteration 1
  ############################################################################################################################
  message(sprintf("%s : %s", Sys.time(), paste0("Running LSI (1 of ",iterations,") on Top Features")))
  j <- 1

  if(!is.null(clusterParams$sampleCells)){
    if(!is.na(clusterParams$sampleCells[j])){
      sampleJ <- clusterParams$sampleCells[j]
    }else if(!is.na(clusterParams$sampleCells[1])){
      sampleJ <- clusterParams$sampleCells[1]
    }else{
      sampleJ <- sampleCellsPre
    }
  }else{
    sampleJ <- sampleCellsPre 
  }

  outLSI <- ArchR:::.LSIPartialMatrix(
    ArrowFiles = ArrowFiles, 
    featureDF = topFeatures,
    cellNames = cellNames, 
    cellDepth = cellDepth,
    useMatrix = useMatrix,
    sampleNames = getCellColData(ArchRProj)$Sample, 
    LSIMethod = LSIMethod, 
    scaleTo = scaleTo,
    dimsToUse = dimsToUse, 
    binarize = binarize, 
    outlierQuantiles = outlierQuantiles,
    sampleCells = if(j != iterations) sampleCellsPre else sampleCellsFinal,
    projectAll = j == iterations | projectCellsPre | sampleJ > sampleCellsPre,
    threads = threads,
    useIndex = FALSE,
    seed = seed,
    tstart = tstart,
    verbose = verbose,
    logFile = logFile
  )
  outLSI$scaleDims <- scaleDims
  outLSI$useMatrix <- useMatrix
  outLSI$tileSize <- tileSize
  gc()

  if(iterations == 1){
    message(sprintf("%s : %s", Sys.time(), "Finished Running IterativeLSI"))
    ArchRProj@reducedDims[[name]] <- outLSI
    return(ArchRProj)
  }

  #########################
  # Identify LSI Clusters
  #########################
  clusterDF <- ArchR:::.LSICluster(
    outLSI = outLSI,
    filterBias = filterBias,
    cellNames = cellNames,
    cellDepth = cellDepth,
    dimsToUse = dimsToUse,
    scaleDims = scaleDims,
    corCutOff = corCutOff,
    clusterParams = clusterParams,
    j = j,
    verbose = verbose,
    tstart = tstart,
    logFile = logFile
  )
  
  if(histology) {message(sprintf("refining cluster using histology"))}
  else {message(sprintf("refining cluster using spot position"))}
                                
  refine_j=1
  while(refine_j <= refine_iter) {
      pred <- cbind(data.frame(clusterDF, row.names = clusterDF$cellNames), data.frame(getCellColData(ArchRProj)[,c("Sample", "array_col", "array_row")]))
      colnames(pred) <- c("cellNames", "clusters", "Sample", "array_row", "array_col")
      if(refine == "default"){
        clusterDF$clusters <- refine(pred = pred, histology = histology)
      } else {
        clusterDF$clusters <- refine(pred = pred, histology = histology)
        pred <- cbind(data.frame(clusterDF, row.names = clusterDF$cellNames), data.frame(getCellColData(ArchRProj)[,c("Sample", "array_col", "array_row")]))
        colnames(pred) <- c("cellNames", "clusters", "Sample", "array_row", "array_col")
        clusterDF$clusters <- refine(pred = pred, histology = histology)
      }
      refine_j=refine_j+1
  }
                                
  clusters <- clusterDF$clusters
  nClust <- length(unique(clusters))
  message(sprintf("Identified %s Clusters", nClust))

  #########################
  # Save LSI Iteration
  #########################
  if(saveIterations){
    message(sprintf("%s : %s", Sys.time(), "Saving LSI Iteration"))
    ArchR:::.saveIteration(outLSI=outLSI, clusters=clusters, scaleDims = scaleDims, 
      dimsToUse = dimsToUse, corCutOff = corCutOff, outDir = outDir,
      nPlot=nPlot, UMAPParams=UMAPParams, ArchRProj=ArchRProj, j = j, threads = threads, logFile = logFile)
  }

  ############################################################################################################################
  # LSI Iteration 2+
  ############################################################################################################################
  variableFeatures <- topFeatures

  while(j < iterations){

    #Jth iteration
    j <- j + 1

    #########################
    # Identify Features for LSI Iteration
    #########################
    variableFeatures <- ArchR:::.identifyVarFeatures(
      outLSI = outLSI,
      clusterDF = clusterDF,
      ArrowFiles = ArrowFiles,
      useMatrix = useMatrix,
      prevFeatures = variableFeatures,
      scaleTo = scaleTo,
      totalAcc = totalAcc,
      totalFeatures = totalFeatures,
      firstSelection = firstSelection,
      selectionMethod = selectionMethod,
      varFeatures = varFeatures,
      tstart = tstart,
      threads = threads,
      verbose = verbose,
      logFile = logFile
    )

    #########################
    # LSI
    #########################
    message(sprintf("Running LSI (%s of %s) on Variable Features", j, iterations))
    if(!is.null(clusterParams$sampleCells)){
      if(!is.na(clusterParams$sampleCells[j])){
        sampleJ <- clusterParams$sampleCells[j]
      }else if(!is.na(clusterParams$sampleCells[1])){
        sampleJ <- clusterParams$sampleCells[1]
      }else{
        sampleJ <- sampleCellsPre
      }
    }else{
      sampleJ <- sampleCellsPre 
    }

    #Compute Partial Matrix LSI
    outLSI <- ArchR:::.LSIPartialMatrix(
      ArrowFiles = ArrowFiles, 
      featureDF = variableFeatures,
      useMatrix = useMatrix,
      cellNames = cellNames, 
      cellDepth = cellDepth,
      sampleNames = getCellColData(ArchRProj)$Sample, 
      LSIMethod = LSIMethod, 
      scaleTo = scaleTo, 
      dimsToUse = dimsToUse,
      binarize = binarize,
      outlierQuantiles = outlierQuantiles, 
      sampleCells = if(j != iterations) sampleCellsPre else sampleCellsFinal,
      projectAll = j == iterations | projectCellsPre | sampleJ > sampleCellsPre,
      threads = threads,
      useIndex = FALSE,
      seed = seed,
      tstart = tstart,
      verbose = verbose,
      logFile = logFile
    )
    outLSI$scaleDims <- scaleDims
    outLSI$useMatrix <- useMatrix
    outLSI$tileSize <- tileSize

    if(j != iterations){

      #########################
      # Identify LSI Clusters
      #########################
      clusterDF <- ArchR:::.LSICluster(
        outLSI = outLSI,
        dimsToUse = dimsToUse,
        scaleDims = scaleDims,
        corCutOff = corCutOff,
        filterBias = filterBias,
        cellNames = cellNames,
        cellDepth = cellDepth,
        j = j,
        clusterParams = clusterParams,
        verbose = verbose,
        tstart = tstart,
        logFile = logFile
      )
        
      
      if(histology) {message(sprintf("refining cluster using histology"))}
      else {message(sprintf("refining cluster using spot position"))}
      
      refine_j=1
      while(refine_j <= refine_iter) {
          pred <- cbind(data.frame(clusterDF, row.names = clusterDF$cellNames), data.frame(getCellColData(ArchRProj)[,c("Sample", "array_col", "array_row")]))
          colnames(pred) <- c("cellNames", "clusters", "Sample", "array_row", "array_col")
          if(refine == "default"){
            clusterDF$clusters <- refine(pred = pred, histology = histology)
          } else {
            clusterDF$clusters <- refine(pred = pred, histology = histology)
            pred <- cbind(data.frame(clusterDF, row.names = clusterDF$cellNames), data.frame(getCellColData(ArchRProj)[,c("Sample", "array_col", "array_row")]))
            colnames(pred) <- c("cellNames", "clusters", "Sample", "array_row", "array_col")
            clusterDF$clusters <- refine(pred = pred, histology = histology)
          }
          refine_j=refine_j+1
      }
      clusters <- clusterDF$clusters
      nClust <- length(unique(clusters))
      message(sprintf("Identified %s Clusters", nClust))

      #########################
      # Save LSI Iteration
      #########################
      if(saveIterations){
        message(sprintf("Saving LSI Iteration"))
        ArchR:::.saveIteration(outLSI=outLSI, clusters=clusters, scaleDims = scaleDims, 
          dimsToUse = dimsToUse, corCutOff = corCutOff, outDir = outDir,
          nPlot=nPlot, UMAPParams=UMAPParams, ArchRProj=ArchRProj, j = j, threads = threads, logFile = logFile)
      }

    }

    gc()

  }

  #Organize Output
  message(sprintf("%s : %s", Sys.time(), "Finished Running IterativeLSI"))
  ArchRProj@reducedDims[[name]] <- outLSI

  return(ArchRProj)

}
                                
refine <- function(pred, num_nbs = 4, histology = FALSE){
    res <- map(unique(pred$Sample), ~{
        sample_id = pred[pred$Sample == .x,]$cellNames
        label = pred[pred$Sample == .x,]$clusters
        x_array = pred[pred$Sample == .x,]$array_row
        y_array = pred[pred$Sample == .x,]$array_col   
        adj_2d <- calculate_adj_matrix(x=x_array,y=y_array)
        refined_pred = refine_modify(sample_id=sample_id, pred=label, dis=adj_2d, shape='square')
        return(refined_pred)
    }) %>% Reduce(c,.)
    return(res)
}

                                
IterativeLSI_Clustering <- function(
  input = NULL, 
  reducedDims = "IterativeLSI",
  name = "Clusters",
  sampleCells = NULL,
  seed = 1, 
  method = "Seurat",
  dimsToUse = NULL,
  scaleDims = NULL, 
  corCutOff = 0.75,
  knnAssign = 10, 
  nOutlier = 5, 
  maxClusters = 25,
  testBias = TRUE,
  filterBias = FALSE,
  biasClusters = 0.01,
  biasCol = "nFrags",
  biasVals = NULL,
  biasQuantiles = c(0.05, 0.95),
  biasEnrich = 10,
  biasProportion = 0.5,
  biasPval = 0.05,
  nPerm = 500,
  prefix = "C",
  ArchRProj = NULL,
  verbose = TRUE,
  tstart = NULL,
  force = FALSE,
  logFile = createLogFile("IterativeLSI_Clustering"),
  ...
  ){

  if(is(ArchRProj, "ArchRProject")){
    message("When running addClusters 'input' param should be used for 'ArchRProj'. Replacing 'input' param with user 'ArchRPRoj'...")
    input <- ArchRProj
    rm(ArchRProj)
    gc()
  }
  
  if(is.null(tstart)){
      tstart <- Sys.time()
  }

  if(inherits(input, "ArchRProject")){
      #Check
      input <- addCellColData(
          ArchRProj = input, 
          data = rep(NA, nCells(input)), 
          name = name, 
          cells = getCellNames(input), 
          force = force
      )

      if(reducedDims %ni% names(input@reducedDims)){
          stop("Error reducedDims not available!")
      }

      matDR <- getReducedDims(
          ArchRProj = input, 
          reducedDims = reducedDims, 
          dimsToUse = dimsToUse, 
          corCutOff = corCutOff, 
          scaleDims = scaleDims
      )
  
  }else if(inherits(input, "matrix")){
      matDR <- input
  }else{
      stop("Input an ArchRProject or Cell by Reduced Dims Matrix!")
  }

  #Subset Matrix
  set.seed(seed)
  nr <- nrow(matDR)

  if(!is.null(sampleCells)){
    if(sampleCells < nrow(matDR)){
      estimatingClusters <- 1
      idx <- sample(seq_len(nrow(matDR)), sampleCells)
      matDRAll <- matDR
      matDR <- matDR[idx,,drop=FALSE]
    }else{
      estimatingClusters <- 0
    }
  }else{
    estimatingClusters <- 0
  }

  #################################################################################
  # Decide on which clustering setup to use
  #################################################################################
  if(grepl("seurat",tolower(method))){
  }else if(grepl("scran",tolower(method))){
  }else{
    stop("Clustering Method Not Recognized!")
  }

  clust <- tryCatch({

    if(grepl("seurat",tolower(method))){

      clustParams <- list(...)
      clustParams$verbose <- verbose
      clustParams$tstart <- tstart
      clust <- ArchR:::.clustSeurat(mat = matDR, clustParams = clustParams, logFile = logFile)

    }else if(grepl("scran",tolower(method))){

      clustParams <- list(...)
      clustParams$verbose <- verbose
      clustParams$tstart <- tstart
      clustParams$x <- t(matDR)
      clustParams$d <- ncol(matDR)
      clustParams$k <- ifelse(exists("...$k"), ...$k, 25)
      clust <- .clustScran(clustParams = clustParams, logFile = logFile)

    }

  }, error = function(e){

    errorList <- clustParams

  })

  #################################################################################
  # If estimating clsuters we will assign to nearest neighbor cluster
  #################################################################################
  if(estimatingClusters == 1){
      
      knnAssigni <- as.matrix(ArchR:::.computeKNN(matDR, matDRAll[-idx,,drop=FALSE], knnAssign))
      clustUnique <- unique(clust)
      clustMatch <- match(clust, clustUnique)
      knnAssigni <- matrix(apply(knnAssigni, 2, function(x) clustMatch[x]), ncol = knnAssign)

      clustAssign <- lapply(seq_along(clustUnique), function(x){
          rowSums(knnAssigni == x)
      }) %>% Reduce("cbind", .) %>% apply(., 1, which.max)
      clustOld <- clust
      clust <- rep(NA, nr)
      clust[idx] <- clustOld
      clust[-idx] <- clustUnique[clustAssign]
      matDR <- matDRAll
      remove(matDRAll)
      gc()

  }

  #################################################################################
  # Testing Bias
  #################################################################################
  if(testBias){
    if(inherits(input, "ArchRProject")){
      if(is.null(biasVals)){
        biasDF <- getCellColData(input, select = biasCol)
      }else{
        biasDF <- DataFrame(row.names = rownames(matDR), bias = biasVals)
      }
    }else{
      if(!is.null(biasVals)){
        biasDF <- DataFrame(row.names = rownames(matDR), bias = biasVals)
      }else{
        message("No biasVals for testing bias continuing without bias detection")
        testBias <- FALSE
      }
    }
  }

  if(testBias){
    clust <- tryCatch({
      biasDF$Q <- ArchR:::.getQuantiles(biasDF[,1])
      tabClust <- table(clust)
      tabClustP <- tabClust / sum(tabClust)
      idxTest <- which(tabClustP < biasClusters)
      names(clust) <- rownames(matDR)
      if(length(idxTest) > 0){
        testDF <- lapply(seq_along(idxTest), function(i){
          clustTesti <- names(tabClustP)[idxTest[i]]
          biasQ <- biasDF[names(clust)[which(clust == clustTesti)], 2]
          biasBgd <- matrix(
            sample(
              x = biasDF[names(clust)[which(clust != clustTesti)], 2],
              size = nPerm * length(biasQ),
              replace = if(nPerm * length(biasQ) > nrow(biasDF[names(clust)[which(clust != clustTesti)], ])) TRUE else FALSE
            ), 
            nrow = length(biasQ), 
            ncol = nPerm
          )
          n1 <- colSums(biasBgd >= max(biasQuantiles))
          n2 <- colSums(biasBgd <= min(biasQuantiles))
          pval1 <- max(sum(sum(biasQ >= max(biasQuantiles)) < n1) * 2, 1) / length(n1)
          pval2 <- max(sum(sum(biasQ <= min(biasQuantiles)) < n2) * 2, 1) / length(n2)
          enrich1 <- sum(biasQ >= max(biasQuantiles)) / max(median(n1), 1)
          enrich2 <- sum(biasQ <= min(biasQuantiles)) / max(median(n2), 1)
          per1 <- sum(biasQ >= max(biasQuantiles)) / length(biasQ)
          per2 <- sum(biasQ <= min(biasQuantiles)) / length(biasQ)
          if(enrich1 > enrich2){
            enrichClust <- enrich1
            enrichPval <- min(pval1, 1)
            enrichPer <- per1
          }else{
            enrichClust <- enrich2
            enrichPval <- min(pval2, 1)
            enrichPer <- per2
          }
          DataFrame(Cluster = clustTesti, enrichClust = enrichClust, enrichPval = enrichPval, enrichProportion = enrichPer)
        }) %>% Reduce("rbind", .)

        clustAssign <- testDF[which(testDF$enrichClust > biasEnrich & testDF$enrichProportion > biasProportion & testDF$enrichPval <= biasPval),1]
        if(length(clustAssign) > 0){
          if(filterBias){
            for(i in seq_along(clustAssign)){
              clusti <- clustAssign[i]
              idxi <- which(clust==clusti)
              knni <- ArchR:::.computeKNN(matDR[-idxi,,drop=FALSE], matDR[idxi,,drop=FALSE], knnAssign)
              clustf <- unlist(lapply(seq_len(nrow(knni)), function(x) names(sort(table(clust[-idxi][knni[x,]]),decreasing=TRUE)[1])))
              clust[idxi] <- clustf
            }
          }else{
            message("Biased Clusters : ", appendLF = FALSE)
            for(i in seq_along(clustAssign)){
              message(clustAssign[i], " ", appendLF = FALSE)
            }
            message("")
          }
        }
      }
      clust
    }, error = function(e){

      errorList <- list(
        idxTest = if(exists("testDF", inherits = FALSE)) fragx else "Error with idxTest!",
        biasDF = if(exists("testDF", inherits = FALSE)) fragx else "Error with biasDF!",
        testDF = if(exists("testDF", inherits = FALSE)) fragx else "Error with testDF!",
        clustAssign = if(exists("idf", inherits = FALSE)) fragx else "Error with clustAssign!"
      )

    })

  }
  
  #################################################################################
  # Test if clusters are outliers identified as cells with fewer than nOutlier
  #################################################################################
  tabClust <- table(clust)
  clustAssign <- which(tabClust < nOutlier)
  if(length(clustAssign) > 0){
      for(i in seq_along(clustAssign)){
          clusti <- names(clustAssign[i])
          idxi <- which(clust==clusti)
          knni <- ArchR:::.computeKNN(matDR[-idxi,], matDR[idxi,], knnAssign)
          clustf <- unlist(lapply(seq_len(nrow(knni)), function(x) names(sort(table(clust[-idxi][knni[x,]]),decreasing=TRUE)[1])))
          clust[idxi] <- clustf
      }
  }

  #################################################################################
  # Merging if more than maxClusters
  #################################################################################
  if(!is.null(maxClusters)){
    if(length(unique(clust)) > maxClusters){
      meanDR <- t(ArchR:::.groupMeans(t(matDR), clust))
      hc <- hclust(dist(as.matrix(meanDR)))
      ct <- cutree(hc, maxClusters)
      clust <- mapLabels(
        labels = clust, 
        oldLabels = names(ct), 
        newLabels = paste0(prefix, ct)
      )
    }
  }
  
  #################################################################################
  # Renaming Clusters based on Proximity in Reduced Dimensions
  #################################################################################
  
  if(length(unique(clust)) > 1){

      meanDR <- t(ArchR:::.groupMeans(t(matDR), clust))
      hc <- hclust(dist(as.matrix(meanDR)))
      out <- mapLabels(
        labels = clust, 
        oldLabels = hc$labels[hc$order], 
        newLabels = paste0(prefix, seq_along(hc$labels))
      )

  }else{

      out <- rep(paste0(prefix, "1"), length(clust))

  }

  if(inherits(input, "ArchRProject")){
    input <- addCellColData(
            input, 
            data = out, 
            name = name, 
            cells = rownames(matDR),
            force = TRUE
        )
    pred <- getCellColData(input)[,c(name,"Sample", "array_col", "array_row")]
    colnames(pred) <- c("clusters", "Sample", "array_row", "array_col")
    pred$cellNames <- rownames(getCellColData(input))
    input@cellColData[,name] <- refine(pred = pred)
    return(input)

  }else if(!inherits(input, "ArchRProject")){
    return(out)
  }

}
