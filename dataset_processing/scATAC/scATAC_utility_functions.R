# This script is in parts taken from 
# https://github.com/GreenleafLab/SpearATAC_MS_2021/blob/main/AnalyzeSpearATAC/Scripts/SpearATAC-Functions.R
# which was published under the MIT licence:
 
# MIT License

# Copyright (c) 2020 GreenleafLab

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.



library(here)

Spear_ATAC_filenames_creator <- function(n_replicates, filename_base, file_type, cell_line){
    replicates <- 1:n_replicates
    type2suffix = c(
        "fragments" = ".fragments.tsv.gz",
        "barcodes" =  ".singlecell.csv",
        "sgRNA_assignment" = ".sgRNA.rds"
    )

    files = c()
    labels = c()
    for (replicate in replicates) {
        file = paste0(filename_base, replicate, type2suffix[file_type], collapse = ", ")
        label = paste0(cell_line, "_", replicate)
        files = c(files, file)
        labels = c(labels, label)
    }
    names(files) <- labels
    return(files)
}


getFileNames <- function(dataset, cell_line, file_type){
    data_dir = paste0(here(), "/scATAC/data/")
    if(dataset == 'Spear_ATAC') {
        data_dir = paste0(data_dir, "Spear_ATAC/")
        if(cell_line == 'K562') {
            filename_base = paste0(data_dir, "K562-LargeScreen-R")
            n_replicates = 6
        }
        if(cell_line == "MCF7") {
            filename_base = paste0(data_dir, "MCF7-LargeScreen-R")
            n_replicates = 4
        }
        if(cell_line == "GM12878") {
            filename_base = paste0(data_dir, "GM-LargeScreen-R")
            n_replicates = 4
        }
        filenames = Spear_ATAC_filenames_creator(n_replicates, filename_base, file_type, cell_line)
    } 
    else if(dataset ==  "CRISPR_sciATAC") {
        data_dir = paste0(data_dir, "CRISPR_sciATAC/")
        if(cell_line == "K562_1") {
          if(file_type == "fragments"){
            filenames = c( paste0(data_dir, "GSM4887677_screen1_snATAC.bed.gz"))
            names(filenames) = "K562"
          }
          else if(file_type == "sgRNA_assignment"){
            filenames = c(paste0(data_dir, "GSM4887678_screen1_snATACguide.IDs.mat.txt.gz"))
            names(filenames) = "K562"
          }
        }
        if(cell_line == "K562_2") {
          if(file_type == "fragments"){
            filenames = c(paste0(data_dir, "GSM4887679_screen2_snATAC.bed.gz"))
            names(filenames) = "K562"
          }
          else if(file_type == "sgRNA_assignment"){
            filenames = c(paste0(data_dir, "GSM4887680_screen2_snATACguide.IDs.mat.txt.gz"))
            names(filenames) = "K562"
          }
        }
    }
    else if(dataset ==  "ASAP_seq") {
        data_dir = paste0(data_dir, "ASAP_seq/")
        if(cell_line == "CD4+_T_cells") {
          if(file_type == "fragments"){
            filenames = c( paste0(data_dir, "GSM4732137_Perturb_CD4_stim_fragments.tsv.gz"))
            names(filenames) = "CD4+_T_cells"
          }
          else if(file_type == "sgRNA_assignment"){
            filenames = c(paste0(data_dir, "HTO_res_filtered.txt"))
            names(filenames) = "CD4+_T_cells"
          }
        }
    }
    return(filenames)
}


#Create an sgRNA assignment matrix
#Iterate over each sgRNA aligned file
assign_sgRNA <- function(sgRNAFiles, proj, dataset){
  sgAssign <- lapply(seq_along(sgRNAFiles), function(x){
    
    message(x)
    if(dataset == "Spear_ATAC"){
      #Read in sgRNA Data Frame
      sgDF <- readRDS(sgRNAFiles[x])   
      #Create sparseMatrix of sgRNA assignments
      sgMat <- createSpMat(sgDF[,1], sgDF[,2])      
      
      #Create Column Names that match EXACTLY with those in the ArchR Project
      #Cell barcodes in our case were the reverse complement to thos in the scATAC-seq data
      colnames(sgMat) <- paste0(names(sgRNAFiles)[x],"#", reverseComplement(DNAStringSet(colnames(sgMat))),"-1")
    }
    if(dataset == "CRISPR_sciATAC"){
      sgMat <- read_csv(sgRNAFiles[x])
      sgMat <- DataFrame(sgMat)
      rownames(sgMat) <- sgMat[,1]
      sgMat <- sgMat[,-1]
      sgMat <- DataFrame(t(data.matrix(sgMat)))
      sgMat <- as.data.frame(sgMat)
      colnames(sgMat) <- paste0(names(sgRNAFiles)[x],"#", colnames(sgMat))
    }
    if(dataset == "ASAP_seq"){

      sg_assign = read_tsv(sgRNAFiles[x] )
      sg_assign[,'cells'] = paste0(names(sgRNAFiles)[x], "#", pull(sg_assign, CellBarcode))
      #filter out cells that don't exist in the proj otherwise there will be trouble later
      sg_assign <- sg_assign[(pull(sg_assign, cells)  %in% proj$cellNames),]

      df <- DataFrame(
        cell = pull(sg_assign, cells),
        sgAssign = pull(sg_assign, HTO),
        sgCounts = Inf, #settings those so that all thresholds are met
        sgTotal = Inf,
        sgSpec = Inf
      )
      return(df)
     
    }
    #Check This
    if(sum(colnames(sgMat) %in% proj$cellNames) < 2){
      stop("x=",x,"; Error matching of sgRNA cell barcodes and scATAC-seq cell barcodes was unsuccessful. Please check your input!")
    }

    #Filter those that are in the ArchR Project
    sgMat <- sgMat[,colnames(sgMat) %in% proj$cellNames,drop=FALSE]

    #Compute sgRNA Staistics
    df <- DataFrame(
      cell = colnames(sgMat), #Cell ID
      sgAssign = rownames(sgMat)[apply(sgMat, 2, which.max)], #Maximum sgRNA counts Assignment
      sgCounts =  apply(sgMat, 2, max), #Number of sgRNA counts for max Assignment
      sgTotal = colSums(sgMat), #Number of total sgRNA counts across all Assignments
      sgSpec = apply(sgMat, 2, max) / colSums(sgMat) #Specificity of sgRNA assignment
    )

    #Return this dataframe
    df

  }) %>% Reduce("rbind", .)

  #Make the rownames the cell barcodes
  rownames(sgAssign) <- sgAssign[,1]
  return(sgAssign)
}






#Access ArchR Hidden Functions
fn <- unclass(lsf.str(envir = asNamespace("ArchR"), all = TRUE))
for (i in seq_along(fn)) {
    tryCatch({
        eval(parse(text = paste0(fn[i], "<-ArchR:::", fn[i])))
    }, error = function(x) {
    })
}

#Function for creating a sparse confusion matrix
createSpMat <- function(i,j){
  ui <- unique(i)
  uj <- unique(j)
  m <- Matrix::sparseMatrix(
    i = match(i, ui),
    j = match(j, uj),
    x = rep(1, length(i)),
    dims = c(length(ui), length(uj))
  )
  rownames(m) <- ui
  colnames(m) <- uj
  m
}

#Function For cleaning sgRNA
cleanSgRNA <- function(
  ArchRProj = NULL, 
  groupSg = NULL,
  individualSg = NULL,
  nonTarget = "sgNT", 
  useMatrix = "TileMatrix",
  scaleTo = 10000, 
  excludeChr = "chrM",
  verbose = FALSE,
  k = 20,
  seed = 1,
  tR  = 0.9,
  varFeatures = 10000,
  dimsToUse = 1:30,
  totalFeatures = 500000,
  filterQuantile = 0.999,
  plotDir = "cleanSgRNA",
  threads = getArchRThreads(),
  logFile = createLogFile("cleanSgRNA")
  ){

  dir.create(plotDir)
  time <- make.names(paste0(Sys.time()))

  fn <- unclass(lsf.str(envir = asNamespace("ArchR"), all = TRUE))
  for (i in seq_along(fn)) {
      tryCatch({
          eval(parse(text = paste0(fn[i], "<-ArchR:::", fn[i])))
      }, error = function(x) {
      })
  }
  
  tstart <- Sys.time()

  .startLogging(logFile=logFile)

  #1 Get Info for LSI

  #MatrixFiles
  ArrowFiles <- getSampleColData(ArchRProj)[,"ArrowFiles"]

  #Check if Matrix is supported and check type
  stopifnot(any(tolower(useMatrix) %in% c("tilematrix","peakmatrix")))
  if(tolower(useMatrix) == "tilematrix"){
    useMatrix <- "TileMatrix"
    tileSizes <- lapply(ArrowFiles, function(x){
      h5read(x, "TileMatrix/Info/Params/")$tileSize[1]
    }) %>% unlist
    if(length(unique(tileSizes)) != 1){
      stop("Error not all TileMatrices are the same tileSize!")
    }
    tileSize <- unique(tileSizes)
  }
  if(tolower(useMatrix) == "peakmatrix"){
    useMatrix <- "PeakMatrix"
    tileSize <- NA
  }

  chrToRun <- .availableSeqnames(ArrowFiles, subGroup = useMatrix)
  
  #Compute Row Sums Across All Samples
  .logDiffTime("Computing Total Accessibility Across All Features", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
  if(useMatrix == "TileMatrix"){
    totalAcc <- .getRowSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun, addInfo = FALSE)
    totalAcc$start <- (totalAcc$idx - 1) * tileSize
  }else{
    totalAcc <- .getRowSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun, addInfo = TRUE)
  }
  gc()

  cellDepth <- tryCatch({
      df <- getCellColData(ArchRProj = ArchRProj, select = "nFrags")
      v <- df[,1]
      names(v) <- rownames(df)
      v
    }, error = function(e){
      tryCatch({
        .getColSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun)
      }, error = function(y){
        stop("Could not determine depth from nFrags or colSums!")
      })
    }
  )
  cellDepth <- log10(cellDepth + 1)

  #Filter Chromosomes
  if(length(excludeChr) > 0){
    totalAcc <- totalAcc[BiocGenerics::which(totalAcc$seqnames %bcni% excludeChr), , drop = FALSE]
  }

  #2 Create sgRNA Group Matrix
  .logDiffTime("Computing Top Features", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
  rmTop <- floor((1-filterQuantile) * totalFeatures)
  topIdx <- head(order(totalAcc$rowSums, decreasing=TRUE), totalFeatures + rmTop)[-seq_len(rmTop)]
  topFeatures <- totalAcc[sort(topIdx),]
  
  #Need to determine which cells are being held out
  sgRNA <- getCellColData(ArchRProj, groupSg)
  groupList <- SimpleList(split(rownames(sgRNA), paste0(sgRNA[,1])))

  groupMat <- .getGroupMatrix(
    ArrowFiles = ArrowFiles, 
    featureDF = topFeatures, 
    groupList = groupList, 
    threads = threads, 
    useIndex = FALSE, 
    verbose = verbose, 
    useMatrix = useMatrix, 
    asSparse = FALSE, 
    tstart = tstart
  )
  groupMat <- t(t(groupMat) / colSums(groupMat)) * scaleTo

  #Get Groups for individual SgRNA
  sgRNA_Individual <- getCellColData(ArchRProj, individualSg)
  groupIndividualList <- SimpleList(split(rownames(sgRNA_Individual), paste0(sgRNA_Individual[,1])))
  
  indMat <- .getGroupMatrix(
    ArrowFiles = ArrowFiles, 
    featureDF = topFeatures, 
    groupList = groupIndividualList, 
    threads = threads, 
    useIndex = FALSE, 
    verbose = verbose, 
    useMatrix = useMatrix, 
    asSparse = FALSE, 
    tstart = tstart
  )
  indMat <- t(t(indMat) / colSums(indMat)) * scaleTo

  #Unique Sg
  uniqSg <- unique(sgRNA[,1])
  uniqSg <- uniqSg[!is.na(uniqSg)]
  uniqSg <- uniqSg[uniqSg %ni% nonTarget]

  mapSg <- split(paste0(sgRNA_Individual[,1]), paste0(sgRNA[,1]))
  mapSg2 <- lapply(seq_along(mapSg), function(x) unique(mapSg[[x]]))
  names(mapSg2) <- names(mapSg)


  #########################################################
  #
  # Individual sgRNA
  #
  #########################################################

  assignList <- lapply(seq_along(uniqSg), function(x){

    message(sprintf("Cleaning sgRNA : %s (%s of %s)", uniqSg[x], x, length(uniqSg)))

    gM <- groupMat[, c(uniqSg[x], nonTarget)]
    gM <- edgeR::cpm(gM, log = TRUE, prior.count = 0.5)

    iM <- indMat[, mapSg2[[uniqSg[x]]],drop=FALSE]
    iM <- edgeR::cpm(iM, log = TRUE, prior.count = 0.5)

    diffUp <- order(gM[,1] - gM[,2], decreasing = TRUE)
    diffUpIdx <- which(rowSums(iM - gM[,2] > 0) > 1)
    featuresUp <- head(diffUp[diffUp %in% diffUpIdx], floor(varFeatures/2))

    diffDo <- order(gM[,2] - gM[,1], decreasing = TRUE)
    diffDoIdx <- which(rowSums(gM[,2] - iM > 0) > 1)
    featuresDo <- head(diffDo[diffDo %in% diffDoIdx], floor(varFeatures/2))

    variableFeatures <- topFeatures[sort(c(featuresUp, featuresDo)),]

  #1. Get Partial Matrix With All Cells of Interest
  mat <- .getPartialMatrix(
    ArrowFiles = ArrowFiles,
    featureDF = variableFeatures,
    useMatrix = useMatrix,
        cellNames = c(groupList[[uniqSg[x]]], groupList[[nonTarget]]),
    doSampleCells = FALSE,
    threads = threads,
    verbose = FALSE
  )  

  #2. Create Initial Manifold
  LSI_T_NT <- .computeLSI(
    mat = mat[,c(groupList[[uniqSg[x]]], groupList[[nonTarget]])], 
    LSIMethod = 2, 
    scaleTo = scaleTo,
    nDimensions = max(dimsToUse),
    binarize = TRUE, 
    verbose = FALSE, 
    tstart = tstart
  )

  #3. Compute UMAP Manifold
  umap_T_NT <- uwot::umap(
    X = .rowZscores(LSI_T_NT[[1]]), 
    n_neighbors = 40, 
    min_dist = 0.4, 
    metric = "cosine", 
    ret_model = TRUE
  )

  #Convert To data.frame
  umap_T_NT <- data.frame(umap_T_NT[[1]])
  rownames(umap_T_NT) <- rownames(LSI_T_NT[[1]])
  colnames(umap_T_NT) <- c("UMAP1", "UMAP2")

  #Compute KNN
  knnIdx <- .computeKNN(umap_T_NT, umap_T_NT, k = k)
  rownames(knnIdx) <- rownames(umap_T_NT)
  targetRatio <- rowSums(knnIdx <= length(groupList[[uniqSg[x]]])) / k

  TPT <- intersect(names(which(targetRatio >= tR)), groupList[[uniqSg[x]]])
  TPNT <- intersect(names(which(targetRatio <= 1 - tR)), groupList[[nonTarget]])

  message(sprintf("TP rate = %s", length(TPT) / length(groupList[[uniqSg[x]]])))

  knnIdx2 <- .computeKNN(umap_T_NT[c(TPT,TPNT),], umap_T_NT, k = k)
  rownames(knnIdx2) <- rownames(umap_T_NT)
  targetRatio2 <- rowSums(knnIdx2 <= length(TPT)) / k

  p0 <- ggPoint(umap_T_NT[,1], umap_T_NT[,2], targetRatio, discrete=F,randomize = TRUE, alpha = 1, xlab="UMAP1",ylab="UMAP2", labelAsFactors=FALSE, labelMeans=FALSE) + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

  p1 <- ggPoint(umap_T_NT[,1], umap_T_NT[,2], targetRatio2, discrete=F,randomize = TRUE, alpha = 1, xlab="UMAP1",ylab="UMAP2", labelAsFactors=FALSE, labelMeans=FALSE) + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

  p2 <- ggPoint(umap_T_NT[,1], umap_T_NT[,2], paste0(seq_len(nrow(umap_T_NT)) <= length(groupList[[uniqSg[x]]])), discrete=T,randomize = TRUE, alpha = 1, xlab="UMAP1",ylab="UMAP2", labelAsFactors=FALSE, labelMeans=FALSE) + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
    pdf(file.path(plotDir, paste0(uniqSg[x], ".pdf")), width = 8, height = 4)
    ggAlignPlots(p0, p1, p2, type="h", draw=TRUE)
    dev.off()

    df <- data.frame(cellNames = names(targetRatio2), targetRatio = targetRatio2, sgRNA = sgRNA[names(targetRatio2),])

  SimpleList(targetRatio = df, UMAP = umap_T_NT)

  })

  df2 <- lapply(seq_along(assignList), function(x){
  assignList[[x]][[1]]
  }) %>% Reduce("rbind", .)

  df2 <- df2[order(df2[,2], decreasing=TRUE),]
  df2 <- df2[!duplicated(paste0(df2[,1])), ]

  sgRNA$PurityRatio <- as.numeric(rep(-1, nrow(sgRNA)))
  sgRNA[paste0(df2[,1]),"PurityRatio"] <- df2[,2]

  sgRNA[which(paste0(sgRNA[,1]) == nonTarget),2] <- 1 - sgRNA[which(paste0(sgRNA[,1]) == nonTarget),2]
  colnames(sgRNA) <- c("sgRNA", "PurityRatio")

  ArchRProj <- addCellColData(ArchRProj, data = sgRNA[,"PurityRatio"], name = "PurityRatio", cells = rownames(sgRNA))

  ArchRProj

}



